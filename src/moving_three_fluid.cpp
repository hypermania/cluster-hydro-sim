#include "moving_three_fluid.hpp"
#include <lapacke.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <iomanip>

namespace {
constexpr double beta=2./3., gamma_gas=5./3.;
double cube(double x) {return x*x*x;}
void require(bool condition,const std::string& message) {
  if(!condition) throw std::runtime_error("MovingThreeFluidSim: "+message);
}
void require(bool condition,const char* message) {
  if(!condition) throw std::runtime_error(std::string("MovingThreeFluidSim: ")+message);
}
}

MovingThreeFluidSim::Vec MovingThreeFluidSim::primitive(int i,int f) const {
  return {value(i,f,RHO),value(i,f,VEL),value(i,f,U)};
}
MovingThreeFluidSim::Vec MovingThreeFluidSim::storage(const Vec& x) const {
  return {x[0],x[0]*x[1],x[0]*(x[2]+param.epsilon*x[1]*x[1]/2)};
}
MovingThreeFluidSim::Vec MovingThreeFluidSim::flux(const Vec& x) const {
  return {x[0]*x[1],x[0]*(x[1]*x[1]+beta*x[2]/param.epsilon),
          x[0]*x[1]*(gamma_gas*x[2]+param.epsilon*x[1]*x[1]/2)};
}
MovingThreeFluidSim::Mat MovingThreeFluidSim::storageJacobian(const Vec& x) const {
  Mat m;
  m<<1,0,0, x[1],x[0],0,
    x[2]+param.epsilon*x[1]*x[1]/2,param.epsilon*x[0]*x[1],x[0];
  return m;
}
MovingThreeFluidSim::Mat MovingThreeFluidSim::fluxJacobian(const Vec& x) const {
  const double r=x[0],v=x[1],u=x[2],e=param.epsilon;
  Mat m;
  m<<v,r,0, v*v+beta*u/e,2*r*v,beta*r/e,
    v*(gamma_gas*u+e*v*v/2),r*(gamma_gas*u+1.5*e*v*v),gamma_gas*r*v;
  return m;
}

MovingThreeFluidSim::Face MovingThreeFluidSim::boundaryFlux(const Vec& x) const {
  Face result;
  const double r=x[0],v=x[1],u=x[2];
  const double c=std::sqrt(gamma_gas*beta*u/param.epsilon);
  if(param.reflecting_boundary) {
    // Odd ghost velocity gives zero mass and energy flux. Freeze acoustic
    // impedance in the implicit wall pressure, as for interior viscosity.
    result.flux[1]=beta*r*u/param.epsilon+r*c*v;
    result.left(1,0)=beta*u/param.epsilon+c*v;
    result.left(1,1)=r*c;result.left(1,2)=beta*r/param.epsilon;
    return result;
  }
  if(v-c>=0) {result.flux=flux(x);result.left=fluxJacobian(x);}
  else if(v+3*c>0) {
    const double chi=(3*c+v)/(4*c);
    Vec fan(r*cube(chi),c*chi,u*chi*chi);
    Mat k;
    k<<cube(chi),3*r*chi*chi/(4*c),-3*r*chi*chi*v/(8*c*u),
       0,0.25,3*c/(8*u), 0,u*chi/(2*c),chi*chi-chi*v/(4*c);
    result.flux=flux(fan);result.left=fluxJacobian(fan)*k;
  }
  return result;
}

void MovingThreeFluidSim::initialize(const std::vector<double>& faces_in,
    const std::vector<double>& initial,bool reference_balance) {
  require(initial.size()%12==0&&initial.size()>=36,"invalid state size");
  state=initial;edges=faces_in;
  const int n=zones(),size=12*n;
  require(edges.size()==static_cast<size_t>(n+1)&&edges[0]==0,"invalid faces");
  centres.resize(n);volumes.resize(n);eta.resize(n);
  for(int i=0;i<n;++i) {
    require(std::isfinite(edges[i+1])&&edges[i+1]>edges[i],"unordered faces");
    centres[i]=(edges[i]+edges[i+1])/2;
    volumes[i]=(cube(edges[i+1])-cube(edges[i]))/3;
    eta[i]=(cube(centres[i])-cube(edges[i]))/(3*volumes[i]);
    for(int f=0;f<NF;++f) state[index(i,f,MASS)]=
      (i?value(i-1,f,MASS):0)+volumes[i]*value(i,f,RHO);
  }
  correction.assign(3*n,0);
  fluxes.resize(3*(n+1));
  band.resize(LDAB*size);rhs.resize(size);factor.resize(LDAB*size);
  increment.resize(size);scale.resize(size);next.resize(size);pivots.resize(size);refinement.resize(size);
  totalTime=0;step=0;Deltat=param.Deltat;
  last_change=linear_residual=last_gravity_work=last_heating=energy_ledger_error=0;
  limiting_index=0;
  cumulative_formed_binaries=cumulative_capture_energy=0;
  last_formed_binaries=last_capture_energy=0;
  escaped_mass.fill(0);escaped_energy.fill(0);escaped_heat.fill(0);
  last_mass_outflow.fill(0);last_energy_outflow.fill(0);last_heat_outflow.fill(0);
  validate();
  if(reference_balance) {
    for(int i=0;i<n;++i) {
      double m=0;for(int f=0;f<NF;++f)
        m+=eta[i]*value(i,f,MASS)+(i?(1-eta[i])*value(i-1,f,MASS):0);
      const double g=m/(centres[i]*centres[i]); // q is intentionally excluded.
      const double al=edges[i]*edges[i],ar=edges[i+1]*edges[i+1];
      for(int f=0;f<NF;++f) {
        const Vec x=primitive(i,f);
        require(x[1]==0,"hydrostatic reference has nonzero velocity");
        // Alternating pressure trace uses the cell to the right of a face.
        const Vec left=x;
        const Vec right=i+1<n?primitive(i+1,f):x;
        correction[3*i+f]=((ar*beta*right[0]*right[2]-al*beta*left[0]*left[2])
          /volumes[i]-(ar-al)/volumes[i]*beta*x[0]*x[2]+x[0]*g)/x[0];
      }
    }
  }
}

void MovingThreeFluidSim::validate() const {
  auto positive=[](double x){return std::isfinite(x)&&x>0;};
  require(zones()>=3&&edges.size()==static_cast<size_t>(zones()+1),"uninitialized geometry");
  require(state.size()%12==0&&state.size()==rhs.size()&&band.size()==LDAB*state.size(),
          "incompatible state/workspace sizes");
  require(positive(param.epsilon)&&std::isfinite(param.q)&&param.q>=0&&
    positive(param.thermal_length_over_radius)&&positive(Deltat)&&
    positive(param.max_timestep)&&positive(param.change_tolerance)&&
    positive(param.StopDensity)&&std::isfinite(param.maxTime)&&param.maxTime>=0&&
    param.maxSteps>=0,"invalid parameters");
  for(double x:param.mass) require(positive(x),"invalid particle mass");
  require((param.reflecting_boundary==0||param.reflecting_boundary==1)&&
    (param.heating_dispersion==HEATING_DONOR_DISPERSION||
     param.heating_dispersion==HEATING_RELATIVE_DISPERSION),"invalid boundary/heating closure");
  require(param.binary_formation==BINARY_FORMATION_OFF||
          param.binary_formation==BINARY_FORMATION_POWER_LAW,"unsupported binary formation mode");
  require(std::isfinite(param.capture_coefficient)&&param.capture_coefficient>=0,
          "invalid capture coefficient");
  if(param.binary_formation==BINARY_FORMATION_POWER_LAW)
    require(std::abs(param.mass[FB]/param.mass[FS]-2)<1e-12,"capture requires m_b=2 m_s");
  for(const auto* a:{&param.c1,&param.c4})for(double x:*a)
    require(std::isfinite(x)&&x>=0,"invalid interaction coefficient");
  for(double x:param.c2)require(std::isfinite(x)&&x>=0,"invalid conductivity");
  for(int i=0;i<zones();++i)for(int f=0;f<NF;++f) {
    if(!(positive(value(i,f,RHO))&&positive(value(i,f,U))&&
      positive(value(i,f,MASS))&&std::isfinite(value(i,f,VEL))))
      require(false,"nonpositive/nonfinite state at cell "+std::to_string(i)+" fluid "+std::to_string(f));
    const double expected=(i?value(i-1,f,MASS):0)+volumes[i]*value(i,f,RHO);
    require(std::abs(value(i,f,MASS)-expected)<1e-9*expected,"mass constraint violated");
  }
}

void MovingThreeFluidSim::buildFluxes() {
  const int n=zones();
  for(int f=0;f<NF;++f) {
    fluxes[f]=Face{};
    for(int j=1;j<n;++j) {
      auto& face=fluxes[3*j+f];face=Face{};
      Vec l=primitive(j-1,f),r=primitive(j,f);
      // No frozen primitive reconstruction: advective diffusion vanishes
      // at rest, and the pressure/source quadrature is balanced separately.
      if(!(l[0]>0&&r[0]>0&&l[2]>0&&r[2]>0))
        require(false,"nonpositive primitive state at face "+std::to_string(j));
      // Advective LLF dissipation. Pressure and enthalpy work stay in the
      // implicit acoustic flux, not in the mass/thermal numerical viscosity.
      const double speed=std::max(std::abs(l[1]),std::abs(r[1]));
      // Alternating acoustic traces: left velocity and right pressure.
      // Their divergence/gradient are an adjoint pair rather than two
      // centred operators with an odd/even null space.
      const double density=(l[0]+r[0])/2,energy=(l[2]+r[2])/2;
      const double velocity=l[1],pressure=beta*r[0]*r[2];
      face.flux=Vec(density*velocity,density*velocity*velocity+pressure/param.epsilon,
        density*velocity*(energy+param.epsilon*velocity*velocity/2)+pressure*velocity);
      face.left<<velocity/2,density,0,
        velocity*velocity/2,2*density*velocity,0,
        velocity*(energy+param.epsilon*velocity*velocity/2)/2,
        density*(energy+1.5*param.epsilon*velocity*velocity)+pressure,density*velocity/2;
      face.right<<velocity/2,0,0,
        velocity*velocity/2+beta*r[2]/param.epsilon,0,beta*r[0]/param.epsilon,
        velocity*(energy+param.epsilon*velocity*velocity/2)/2+beta*r[2]*velocity,
        0,density*velocity/2+beta*r[0]*velocity;
      face.flux-=speed*(storage(r)-storage(l))/2;
      face.left+=speed*storageJacobian(l)/2;
      face.right-=speed*storageJacobian(r)/2;
      // Frozen velocity viscosity damps the collocated acoustic mode without
      // diffusing hydrostatic density/entropy at the sound speed. Include its
      // mechanical energy flux; diffusing rho*v at acoustic speed instead
      // produced a trace-fluid odd/even instability on changing backgrounds.
      const double acoustic=std::max(std::abs(l[1])+std::sqrt(gamma_gas*beta*l[2]/param.epsilon),
                                     std::abs(r[1])+std::sqrt(gamma_gas*beta*r[2]/param.epsilon));
      const double viscosity=(acoustic-speed)*(l[0]+r[0])/2;
      face.flux[1]-=viscosity*(r[1]-l[1])/2;
      face.flux[2]-=param.epsilon*viscosity*(r[1]*r[1]-l[1]*l[1])/4;
      face.left(1,1)+=viscosity/2;face.right(1,1)-=viscosity/2;
      face.left(2,1)+=param.epsilon*viscosity*l[1]/2;
      face.right(2,1)-=param.epsilon*viscosity*r[1]/2;
      const double rl=value(j-1,f,RHO),rr=value(j,f,RHO);
      const double k=param.c2[f]*edges[j]*edges[j]*2*rl*rr/(rl+rr)/(centres[j]-centres[j-1]);
      const double zl=std::sqrt(value(j-1,f,U)),zr=std::sqrt(value(j,f,U));
      face.heat=-k*(zr-zl);face.heat_left=k/(2*zl);face.heat_right=-k/(2*zr);
    }
    auto& face=fluxes[3*n+f];face=boundaryFlux(primitive(n-1,f));
    if(param.reflecting_boundary)continue;
    const double h=param.c2[f]*edges[n]*edges[n]*value(n-1,f,RHO)/
      (param.thermal_length_over_radius*edges[n]+edges[n]-centres[n-1]);
    face.heat=h*std::sqrt(value(n-1,f,U));
    face.heat_left=h/(2*std::sqrt(value(n-1,f,U)));
  }
}

void MovingThreeFluidSim::add(int row,int col,double x) {
  if(x==0)return;
  require(row-col<=KL&&col-row<=KU,"coefficient outside optimized band");
  band[KL+KU+row-col+LDAB*col]+=x;
}

void MovingThreeFluidSim::thermalSource(int i,Vec& source,Mat& derivative) const {
  source.setZero();derivative.setZero();
  for(int f=0;f<NF;++f)for(int h=0;h<NF;++h) {
    const double r=value(i,f,RHO),u=value(i,f,U);
    if(h!=f) {
      const double uh=value(i,h,U);
      const double coeff=r*param.c1[3*f+h]*value(i,h,RHO)/
        (param.mass[0]*std::pow(u+uh,1.5));
      source[f]-=coeff*(param.mass[f]*u-param.mass[h]*uh);
      derivative(f,f)-=coeff*param.mass[f];derivative(f,h)+=coeff*param.mass[h];
    }
    if(param.heating_dispersion==HEATING_RELATIVE_DISPERSION) {
      const double sum=u+value(i,h,U),b=r*param.c4[3*f+h]*value(i,h,RHO);
      source[f]+=b/std::sqrt(sum);
      const double slope=-b/(2*std::pow(sum,1.5));
      derivative(f,f)+=slope;derivative(f,h)+=slope;
    } else if((f==FS&&h==FB)||(f==FB)||(f==FD&&h==FB)) {
      const int donor=f==FB?h:f;
      const double ud=value(i,donor,U),b=r*param.c4[3*f+h]*value(i,h,RHO);
      source[f]+=b/std::sqrt(ud);derivative(f,donor)-=b/(2*std::pow(ud,1.5));
    }
  }
}

double MovingThreeFluidSim::captureFrequency(int i) const {
  return param.binary_formation==BINARY_FORMATION_OFF ? 0 :
    param.mass[FB]*param.capture_coefficient*value(i,FS,RHO)/std::pow(value(i,FS,U),0.6);
}

void MovingThreeFluidSim::formationSource(int i,int f,Vec& source,Mat& derivative) const {
  source.setZero();derivative.setZero();
  if(f==FD||param.binary_formation==BINARY_FORMATION_OFF)return;
  const Vec x=primitive(i,FS);
  const double k=captureFrequency(i);
  source=storage(x);derivative=storageJacobian(x);
  if(f==FB) {
    source[2]-=0.5*x[0]*x[2];
    derivative(2,0)-=0.5*x[2];derivative(2,2)-=0.5*x[0];
  }
  const double sign=f==FS?-1:1;
  source*=sign*k;derivative*=sign*k;
}

void MovingThreeFluidSim::assembleCells() {
  const int n=zones();const double dt=Deltat,e=param.epsilon;
  for(int i=0;i<n;++i) {
    const double radius2=centres[i]*centres[i],vol=volumes[i];
    const double al=edges[i]*edges[i],ar=edges[i+1]*edges[i+1],a=(ar-al)/vol;
    double m=0;for(int f=0;f<NF;++f)
      m+=eta[i]*value(i,f,MASS)+(i?(1-eta[i])*value(i-1,f,MASS):0);
    const double gravity=m/radius2-param.q*centres[i];
    Vec heat;Mat heat_derivative;thermalSource(i,heat,heat_derivative);
    for(int f=0;f<NF;++f) {
      const Vec x=primitive(i,f);const double r=x[0],v=x[1],u=x[2];
      // Reference quadrature correction is an acceleration, not a constant
      // force density acting on vanishing gas. Include its mechanical work
      // in energy as well, or expansion spuriously destroys internal energy.
      const double g=gravity-correction[3*i+f];
      Vec source(0,(a*beta*r*u-r*g)/e,-r*v*g);
      source[2]+=heat[f];
      Vec formation;Mat formation_derivative;
      formationSource(i,f,formation,formation_derivative);
      source+=formation;
      for(int row=0;row<3;++row)for(int col=0;col<3;++col)
        add(index(i,f,row),index(i,FS,col),-dt*vol*formation_derivative(row,col));
      Mat deriv=Mat::Zero();
      deriv(1,0)=(a*beta*u-g)/e;deriv(1,2)=a*beta*r/e;
      deriv(2,0)=-v*g;deriv(2,1)=-r*g;
      Mat local=vol*(storageJacobian(x)-dt*deriv);
      for(int h=0;h<NF;++h) {
        add(index(i,f,U),index(i,h,U),-dt*vol*heat_derivative(f,h));
        for(int side=0;side<2;++side) if(i||side==0) {
          const double weight=side?1-eta[i]:eta[i];
          add(index(i,f,VEL),index(i-side,h,MASS),dt*vol*r*weight/(e*radius2));
          add(index(i,f,U),index(i-side,h,MASS),dt*vol*r*v*weight/radius2);
        }
      }
      for(int arow=0;arow<3;++arow) {
        const int row=index(i,f,arow);
        rhs[row]=dt*vol*source[arow];
        for(int k=0;k<3;++k)add(row,index(i,f,k),local(arow,k));
        for(int side=0;side<2;++side) {
          const int j=i+side;const auto& face=fluxes[3*j+f];
          const double sign=side?1:-1,area=edges[j]*edges[j];
          rhs[row]-=dt*sign*area*face.flux[arow];
          for(int k=0;k<3;++k) {
            if(j>0)add(row,index(j-1,f,k),dt*sign*area*face.left(arow,k));
            if(j<n)add(row,index(j,f,k),dt*sign*area*face.right(arow,k));
          }
          if(arow==U) {
            rhs[row]-=dt*sign*face.heat;
            if(j>0)add(row,index(j-1,f,U),dt*sign*face.heat_left);
            if(j<n)add(row,index(j,f,U),dt*sign*face.heat_right);
          }
        }
      }
      const int row=index(i,f,MASS);
      add(row,row,1);add(row,index(i,f,RHO),-vol);
      if(i)add(row,index(i-1,f,MASS),-1);
      rhs[row]=-(value(i,f,MASS)-(i?value(i-1,f,MASS):0)-vol*r);
    }
  }
}

void MovingThreeFluidSim::assembleStep() {
  require(zones()>=3&&state.size()%12==0&&state.size()==rhs.size()&&
          band.size()==LDAB*state.size(),"incompatible state/workspace sizes");
  require(std::isfinite(Deltat)&&Deltat>0&&totalTime+Deltat>totalTime,"invalid timestep");
  std::fill(band.begin(),band.end(),0);std::fill(rhs.begin(),rhs.end(),0);
  buildFluxes();assembleCells();
}

void MovingThreeFluidSim::solveBanded() {
  const int size=static_cast<int>(state.size());
  // Column scaling uses local physical magnitudes, followed by max-row
  // equilibration. It changes neither the equations nor their bandwidth.
  for(int i=0;i<zones();++i)for(int f=0;f<NF;++f)for(int k=0;k<4;++k)
    scale[index(i,f,k)]=k==VEL?std::sqrt(value(i,f,U)/param.epsilon):value(i,f,k);
  factor=band;increment=rhs;std::fill(next.begin(),next.end(),0);
  for(int col=0;col<size;++col)for(int row=std::max(0,col-KU);row<=std::min(size-1,col+KL);++row) {
    double& x=factor[KL+KU+row-col+LDAB*col];x*=scale[col];
    next[row]=std::max(next[row],std::abs(x));
  }
  for(int row=0;row<size;++row) {
    require(std::isfinite(next[row])&&next[row]>0,"invalid matrix row");
    increment[row]/=next[row];
  }
  for(int col=0;col<size;++col)for(int row=std::max(0,col-KU);row<=std::min(size-1,col+KL);++row)
    factor[KL+KU+row-col+LDAB*col]/=next[row];
  const int info=LAPACKE_dgbsv(LAPACK_COL_MAJOR,size,KL,KU,1,factor.data(),LDAB,
                               pivots.data(),increment.data(),size);
  require(info==0,"band solve failed, info="+std::to_string(info));
  for(int k=0;k<size;++k)increment[k]*=scale[k];
  int worst_row=0;
  // Source transfer can leave another component's nearly homogeneous rows
  // inaccurate despite a good normwise solve. Refine against the ORIGINAL
  // equations, reusing the band LU and its existing row/column scales.
  for(int pass=0;pass<4;++pass) {
    linear_residual=0;
    for(int row=0;row<size;++row) {
      long double residual=-rhs[row],norm=std::abs(rhs[row]);
      for(int col=std::max(0,row-KL);col<=std::min(size-1,row+KU);++col) {
        const long double term=static_cast<long double>(band[KL+KU+row-col+LDAB*col])*increment[col];
        residual+=term;norm+=std::abs(term);
      }
      refinement[row]=-residual/next[row];
      const double relative=std::abs(residual)/(norm+1e-300L);
      require(std::isfinite(relative),"nonfinite linear residual");
      if(relative>linear_residual) {linear_residual=relative;worst_row=row;}
    }
    if(linear_residual<1e-8||pass==3)break;
    const int status=LAPACKE_dgbtrs(LAPACK_COL_MAJOR,'N',size,KL,KU,1,factor.data(),LDAB,
                                    pivots.data(),refinement.data(),size);
    require(status==0,"band refinement failed");
    for(int k=0;k<size;++k)increment[k]+=scale[k]*refinement[k];
  }
  if(!(std::isfinite(linear_residual)&&linear_residual<1e-8)) {
    std::ostringstream message;message<<std::setprecision(17)<<"large linear backward error "
      <<linear_residual<<" at row "<<worst_row<<" step "<<step<<" dt "<<Deltat;
    require(false,message.str());
  }
}

double MovingThreeFluidSim::recoverAndBudget() {
  double change=0;
  long double energy_before=0,energy_after=0;
  last_gravity_work=0;last_heating=0;
  last_formed_binaries=last_capture_energy=0;
  for(int i=0;i<zones();++i) {
    double m=0,dm=0;
    for(int f=0;f<NF;++f) {
      m+=eta[i]*value(i,f,MASS)+(i?(1-eta[i])*value(i-1,f,MASS):0);
      dm+=eta[i]*increment[index(i,f,MASS)]+(i?(1-eta[i])*increment[index(i-1,f,MASS)]:0);
    }
    const double gravity=m/(centres[i]*centres[i])-param.q*centres[i];
    const double dg=dm/(centres[i]*centres[i]);
    Vec heat;Mat derivative;thermalSource(i,heat,derivative);
    Vec du;for(int f=0;f<NF;++f)du[f]=increment[index(i,f,U)];
    last_heating+=Deltat*volumes[i]*(heat+derivative*du).sum();
    const double capture=captureFrequency(i);
    if(capture>0) {
      const double r=value(i,FS,RHO),u=value(i,FS,U);
      const double dr=increment[index(i,FS,RHO)],du=increment[index(i,FS,U)];
      const double transfer=Deltat*volumes[i]*capture*(r+dr);
      require(std::isfinite(transfer)&&transfer>=0,"invalid implicit capture transfer");
      last_formed_binaries+=transfer/param.mass[FB];
      last_capture_energy-=0.5*Deltat*volumes[i]*capture*(r*u+u*dr+r*du);
    }
    for(int f=0;f<NF;++f) {
      const double r=value(i,f,RHO),v=value(i,f,VEL),g=gravity-correction[3*i+f];
      last_gravity_work-=Deltat*volumes[i]*(r*v*g+v*g*increment[index(i,f,RHO)]+
        r*g*increment[index(i,f,VEL)]+r*v*dg);
    }
  }
  for(int i=0;i<zones();++i)for(int f=0;f<NF;++f) {
    const Vec x=primitive(i,f);
    const Vec dx(increment[index(i,f,RHO)],increment[index(i,f,VEL)],increment[index(i,f,U)]);
    const Vec w=storage(x)+storageJacobian(x)*dx;
    const double rho=w[0],v=w[1]/rho,u=w[2]/rho-param.epsilon*v*v/2;
    energy_before+=static_cast<long double>(volumes[i])*storage(x)[2];
    energy_after+=static_cast<long double>(volumes[i])*rho*(u+param.epsilon*v*v/2);
    if(!(std::isfinite(rho)&&rho>0&&std::isfinite(v)&&std::isfinite(u)&&u>0))
      require(false,"conservative recovery lost positivity at cell "+std::to_string(i)+" fluid "+std::to_string(f));
    next[index(i,f,RHO)]=rho;next[index(i,f,VEL)]=v;next[index(i,f,U)]=u;
    next[index(i,f,MASS)]=value(i,f,MASS)+increment[index(i,f,MASS)];
    for(int k=0;k<3;++k) {
      const double error=k==RHO?std::abs(rho/x[0]-1):(k==U?std::abs(u/x[2]-1):
        std::abs(v-x[1])*std::sqrt(param.epsilon/x[2]));
      if(error>change) {change=error;limiting_index=index(i,f,k);}
    }
  }
  for(int f=0;f<NF;++f) {
    const int i=zones()-1;const auto& face=fluxes[3*zones()+f];
    const Vec dx(increment[index(i,f,RHO)],increment[index(i,f,VEL)],increment[index(i,f,U)]);
    const Vec out=face.flux+face.left*dx;
    const double area=edges.back()*edges.back();
    last_mass_outflow[f]=Deltat*area*out[0];
    last_energy_outflow[f]=Deltat*area*out[2];
    last_heat_outflow[f]=Deltat*(face.heat+face.heat_left*dx[2]);
    require(last_mass_outflow[f]>=-1e-12*value(i,f,MASS),"linearized boundary caused inflow");
    escaped_mass[f]+=last_mass_outflow[f];escaped_energy[f]+=last_energy_outflow[f];
    escaped_heat[f]+=last_heat_outflow[f];
  }
  long double residual=energy_after-energy_before-last_gravity_work-last_heating-last_capture_energy;
  for(int f=0;f<NF;++f)residual+=last_energy_outflow[f]+last_heat_outflow[f];
  energy_ledger_error=std::abs(residual)/std::max(energy_before,1e-300L);
  require(std::isfinite(energy_ledger_error)&&energy_ledger_error<1e-8,"energy/work ledger failed");
  cumulative_formed_binaries+=last_formed_binaries;
  cumulative_capture_energy+=last_capture_energy;
  state.swap(next);last_change=change;return change;
}

void MovingThreeFluidSim::selectNextTimestep(double change) {
  const double factor=change>0?std::min(1.25,param.change_tolerance/change):1.25;
  Deltat=std::min(param.max_timestep,Deltat*factor);
  require(std::isfinite(Deltat)&&Deltat>0,"invalid next timestep");
}
void MovingThreeFluidSim::advanceAcceptedStep() {
  assembleStep();solveBanded();const double change=recoverAndBudget();
  totalTime+=Deltat;++step;validate();selectNextTimestep(change);
}
double MovingThreeFluidSim::centralDensity() const {
  return std::max({value(0,FS,RHO),value(0,FB,RHO),value(0,FD,RHO)});
}
bool MovingThreeFluidSim::stopCondition() const {
  return totalTime>=param.maxTime||step>=param.maxSteps||centralDensity()>=param.StopDensity;
}
