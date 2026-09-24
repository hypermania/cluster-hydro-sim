#include "../src/moving_initialization.hpp"
#include "../src/moving_comparison.hpp"
#include "../src/moving_statler_observer.hpp"
#include "../src/statler_reproduction.hpp"
#include "../src/moving_statler_initialization.hpp"
#include "../src/heggie_reproduction.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <cstring>
#include <chrono>

using S=MovingThreeFluidSim;
void check(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
S setup(int n=32) {
  ThreeFluidSim hydro;hydro.initSolver(n);hydro.initPlummer(1,0.1,0.2,1,1);
  S s;s.param.epsilon=0.03;s.param.Deltat=1e-5;
  initializeMovingFromHydrostatic(s,hydro);return s;
}
Eigen::MatrixXd dense(const S& s) {
  int n=s.state.size();Eigen::MatrixXd a=Eigen::MatrixXd::Zero(n,n);
  for(int j=0;j<n;++j)for(int i=std::max(0,j-S::KU);i<=std::min(n-1,j+S::KL);++i)
    a(i,j)=s.bandMatrix()[S::KL+S::KU+i-j+S::LDAB*j];
  return a;
}
void bandSolve() {
  auto s=setup(12);s.param.c2={0.1,0.03,0.02};s.param.c1.fill(0.04);s.param.c4.fill(0.005);
  s.param.q=0.002;s.assembleStep();
  const auto a=dense(s);const auto rhs=Eigen::Map<const Eigen::VectorXd>(s.rightHandSide().data(),s.state.size()).eval();
  const Eigen::VectorXd expected=a.partialPivLu().solve(rhs);
  s.advanceAcceptedStep();const auto got=Eigen::Map<const Eigen::VectorXd>(s.lastIncrement().data(),s.state.size());
  check((a*got-rhs).norm()/(rhs.norm()+1e-300)<1e-10,"band residual");
  check((got-expected).norm()/(expected.norm()+1e-300)<1e-8,"band/dense disagreement");
}
void massLedger() {
  auto s=setup();std::array<double,3> mass;
  for(int f=0;f<3;++f)mass[f]=s.value(s.zones()-1,f,S::MASS);
  s.param.maxSteps=40;s.param.maxTime=1;s.param.c2={0.1,0.03,0.02};
  int calls=0;auto observer=[&](const S& x) {
    check(x.energy_ledger_error<1e-12,"energy/work ledger");
    ++calls;for(int f=0;f<3;++f)check(std::abs(x.value(x.zones()-1,f,S::MASS)+x.escaped_mass[f]-mass[f])<1e-12*mass[f],"mass ledger");
  };s.evolve(observer);check(calls==41,"observer calls");
  check(s.escaped_mass[0]>0&&s.escaped_heat[0]>0,"missing open fluxes");
}
void captureLedger() {
  auto zero=setup(12),off_zero=zero;
  zero.param.binary_formation=BINARY_FORMATION_POWER_LAW;
  zero.advanceAcceptedStep();off_zero.advanceAcceptedStep();
  check(zero.state==off_zero.state,"zero capture changes evolution");
  auto s=setup(24);s.param.binary_formation=BINARY_FORMATION_POWER_LAW;
  s.param.capture_coefficient=1e4;s.param.maxSteps=40;
  std::array<double,3> initial;
  for(int f=0;f<NF;++f)initial[f]=s.value(s.zones()-1,f,S::MASS);
  auto observer=[&](const S& x) {
    check(x.energy_ledger_error<1e-12,"capture energy ledger");
    for(int f=0;f<NF;++f) {
      const double transfer=(f==FS?-1:f==FB?1:0)*x.param.mass[FB]*x.cumulative_formed_binaries;
      const double error=x.value(x.zones()-1,f,S::MASS)+x.escaped_mass[f]-initial[f]-transfer;
      check(std::abs(error)<1e-12*initial[f],"capture species mass ledger");
    }
  };
  s.evolve(observer);
  check(s.cumulative_formed_binaries>0&&s.cumulative_capture_energy<0,"missing capture ledgers");

  // Matrix-only audit at nonzero velocity: formation must transfer both
  // momentum and bulk energy, not just rho/U. Compare to an OFF twin.
  auto off=setup(12);
  for(int i=0;i<off.zones();++i)off.state[S::index(i,FS,S::VEL)]=0.17;
  auto on=off;on.param.binary_formation=BINARY_FORMATION_POWER_LAW;
  on.param.capture_coefficient=1e4;off.assembleStep();on.assembleStep();
  const auto da=(dense(on)-dense(off)).eval();
  for(int i=0;i<on.zones();++i) {
    const double r=on.value(i,FS,S::RHO),u=on.value(i,FS,S::U),v=on.value(i,FS,S::VEL);
    const double k=on.param.mass[FB]*on.param.capture_coefficient*r/std::pow(u,0.6);
    const double dtvol=on.Deltat*on.volume(i);
    for(int row=0;row<3;++row) {
      const double combined=on.rightHandSide()[S::index(i,FS,row)]+on.rightHandSide()[S::index(i,FB,row)]
        -off.rightHandSide()[S::index(i,FS,row)]-off.rightHandSide()[S::index(i,FB,row)];
      const double expected=row==S::U?-0.5*dtvol*k*r*u:0;
      check(std::abs(combined-expected)<1e-13*(1+std::abs(expected)),"formation RHS sum");
    }
    check(std::abs(da(S::index(i,FB,S::VEL),S::index(i,FS,S::RHO))+dtvol*k*v)<1e-14,
      "formation momentum derivative");
    check(std::abs(da(S::index(i,FB,S::U),S::index(i,FS,S::U))+0.5*dtvol*k*r)<1e-14,
      "formation birth energy derivative");
  }
  on.advanceAcceptedStep();
  check(on.energy_ledger_error<1e-12,"moving capture kinetic-energy ledger");
}

void movingStatlerSmoke() {
  ThreeFluidSim h;h.param=statlerParameters();h.param.N=32;
  const StatlerInitParam init;initializeCaptureCluster(h,init);
  S s;s.param.Deltat=1e-9;s.param.maxSteps=20;initializeMovingFromHydrostatic(s,h);
  const std::string dir="output/check_moving_statler/";
  std::filesystem::create_directories(dir);
  MovingStatlerObserver observer(statlerObserverParameters(init),dir);
  s.evolve(observer);observer.save(dir);
  check(s.step==20&&s.cumulative_formed_binaries>0,"formation not reached through evolve");
  check(observer.diagnostic.history.time_trh.size()==21,"moving Statler observer calls");
  check(observer.diagnostic.snapshots.snapshots[0].radius[0]==s.radius(0),"moving snapshot cell centres");
  MovingStatlerGrid grid;grid.zones=80;
  S matched;initializeMovingCaptureCluster(matched,init,grid);
  check(matched.faces()[1]==grid.first_face&&std::abs(matched.faces().back()/grid.last_face-1)<1e-14,
        "Statler radial domain");
  double mass=0;for(int f=0;f<NF;++f)mass+=matched.value(matched.zones()-1,f,S::MASS);
  check(std::abs(mass/matched.param.mass[FS]/init.stellar_number-1)<1e-14,"matched stellar number");
  check(std::abs(matched.value(0,FS,S::U)*12-1)<1e-10,"analytic Plummer energy units");
  matched.param.maxSteps=10;auto no_output=[](const S&){};matched.evolve(no_output);
  check(matched.cumulative_formed_binaries>0,"matched-domain capture");
}
void equilibriumAndTide() {
  auto s=setup(80);s.assembleStep();
  for(int i=0;i<s.zones()-1;++i)for(int f=0;f<3;++f) {
    check(std::abs(s.rightHandSide()[S::index(i,f,S::RHO)])<1e-20,"equilibrium mass residual");
    check(std::abs(s.rightHandSide()[S::index(i,f,S::U)])<1e-20,"equilibrium energy residual");
    check(std::abs(s.rightHandSide()[S::index(i,f,S::VEL)])<1e-14,"equilibrium momentum residual");
  }
  const auto original=s.rightHandSide();s.param.q=2;s.assembleStep();
  for(int i=0;i<s.zones();++i)for(int f=0;f<3;++f) {
    const double expected=s.Deltat*s.volume(i)*s.value(i,f,S::RHO)*s.param.q*s.radius(i)/s.param.epsilon;
    const double got=s.rightHandSide()[S::index(i,f,S::VEL)]-original[S::index(i,f,S::VEL)];
    check(std::abs(got-expected)<1e-12*expected,"tidal source or force-balance guard");
  }
}
void reflectingAndRelativeHeating() {
  auto s=setup(24);s.param.reflecting_boundary=1;s.param.c2.fill(.1);
  s.param.maxSteps=30;
  std::array<double,3> initial{};
  for(int f=0;f<NF;++f)initial[f]=s.value(s.zones()-1,f,S::MASS);
  auto observer=[&](const S& x) {
    for(int f=0;f<NF;++f) {
      check(x.escaped_mass[f]==0&&x.escaped_energy[f]==0&&x.escaped_heat[f]==0,"reflecting wall leaked");
      check(std::abs(x.value(x.zones()-1,f,S::MASS)/initial[f]-1)<1e-12,"reflecting mass conservation");
    }
    check(x.energy_ledger_error<1e-12,"reflecting energy ledger");
  };s.evolve(observer);
  auto off=setup(12),on=off;
  on.param.heating_dispersion=HEATING_RELATIVE_DISPERSION;
  on.param.c4[FS*3+FB]=.6;off.assembleStep();on.assembleStep();
  const auto difference=(dense(on)-dense(off)).eval();
  for(int i=0;i<4;++i) {
    const double sum=on.value(i,FS,S::U)+on.value(i,FB,S::U);
    const double b=.6*on.value(i,FS,S::RHO)*on.value(i,FB,S::RHO);
    const double source=on.Deltat*on.volume(i)*b/std::sqrt(sum);
    const int row=S::index(i,FS,S::U);
    check(std::abs((on.rightHandSide()[row]-off.rightHandSide()[row])/source-1)<1e-12,
          "relative dispersion heating RHS");
    for(int f:{FS,FB})check(std::abs(difference(row,S::index(i,f,S::U))/(source/(2*sum))-1)<1e-10,
          "relative dispersion cross-temperature derivative");
  }
  on=off;on.param.heating_dispersion=HEATING_RELATIVE_DISPERSION;
  on.param.c4[FB*3+FB]=.6;on.assembleStep();
  const auto diagonal=(dense(on)-dense(off)).eval();
  for(int i=0;i<4;++i) {
    const int row=S::index(i,FB,S::U);
    const double u=on.value(i,FB,S::U),r=on.value(i,FB,S::RHO);
    const double source=on.Deltat*on.volume(i)*.6*r*r/std::sqrt(2*u);
    check(std::abs(diagonal(row,row)/(source/(2*u))-1)<1e-10,
          "same-component heating must include both temperature derivatives");
  }
}

void heggieSmoke() {
  for(int model=0;model<4;++model) {
    HeggieInitParam initial;initial.zones=32;initial.model=model;
    S s;s.param.maxSteps=20;const auto units=initializeHeggie(s,initial);
    check(s.param.binary_formation==BINARY_FORMATION_OFF&&s.param.reflecting_boundary==1,
          "Heggie physics selection");
    if(model==1)check(std::abs(s.value(0,FB,S::RHO)/s.value(0,FS,S::RHO)-.01)<1e-15,
                     "Heggie segregation fraction");
    if(model>=2)check(s.param.c4[FS*3+FB]==2*s.param.c4[FB*3+FS],"Heggie heating partition");
    check(units.time_unit_over_trh>2&&units.time_unit_over_trh<3,"Heggie time conversion");
    const std::string dir="output/check_heggie_"+std::to_string(model)+"/";
    std::filesystem::create_directories(dir);HeggieObserver observer(units,dir);
    s.evolve(observer);observer.finish(dir,s,"SMOKE");
    check(std::filesystem::file_size(dir+"history.dat")==21*13*sizeof(double),"Heggie binary history layout");
  }
}
void splitSymmetry() {
  auto s=setup(48);for(int i=0;i<s.zones();++i) {
    s.state[S::index(i,FD,S::RHO)]=s.value(i,FS,S::RHO);
    s.state[S::index(i,FD,S::U)]=s.value(i,FS,S::U);
  }
  auto faces=s.faces();auto state=s.state;s.initialize(faces,state,true);
  s.param.mass[FD]=s.param.mass[FS];s.param.c2[FD]=s.param.c2[FS]=0.17;
  s.param.c1[FS*3+FD]=s.param.c1[FD*3+FS]=0.325;
  s.param.maxSteps=25;auto obs=[](const S&){};s.evolve(obs);
  for(int i=0;i<s.zones();++i)for(int k:{S::RHO,S::U,S::MASS})
    check(std::abs(s.value(i,FS,k)/s.value(i,FD,k)-1)<1e-10,"equal-fluid symmetry");
}
void failEarly() {
  auto s=setup();s.param.epsilon=0;
  bool failed=false;try{s.validate();}catch(const std::exception&){failed=true;}
  check(failed,"epsilon validation");
  s=setup();s.Deltat=0;failed=false;
  try{s.advanceAcceptedStep();}catch(const std::exception&){failed=true;}
  check(failed&&s.step==0,"invalid timestep accepted");
  s=setup();s.state.push_back(1);failed=false;
  try{s.advanceAcceptedStep();}catch(const std::exception&){failed=true;}
  check(failed,"resized state reached band solve");
  for(long long mode:{BINARY_FORMATION_MODE_1,BINARY_FORMATION_MODE_2,99LL}) {
    s=setup();s.param.binary_formation=mode;failed=false;
    try{s.validate();}catch(const std::exception&){failed=true;}
    check(failed,"unsupported capture mode accepted");
  }
  s=setup();s.param.capture_coefficient=-1;failed=false;
  try{s.validate();}catch(const std::exception&){failed=true;}
  check(failed,"negative capture coefficient accepted");
}
void dilutedReference() {
  auto s=setup();
  for(int i=s.zones()-4;i<s.zones();++i)for(int f=0;f<3;++f) {
    s.state[S::index(i,f,S::RHO)]*=1e-4;
    s.state[S::index(i,f,S::U)]*=0.01;
  }
  // This used to fail while constructing a face, before any solve: the
  // old equilibrium offset subtracted more density than remained locally.
  s.assembleStep();
  for(double x:s.bandMatrix())check(std::isfinite(x),"diluted reference coefficient");
}
void tidalEvolution() {
  auto s=setup();s.param.q=2;s.param.maxSteps=10;
  auto obs=[](const S&){};s.evolve(obs);
  check(s.value(0,FS,S::VEL)>0&&s.escaped_mass[0]>0,"outward q evolution");
}
void thermalAgreement() {
  ThreeFluidSim h;h.initSolver(32);h.initPlummer(1,0.1,0.2,1,1);
  h.param.md=0.1*h.param.ms;h.param.c1.fill(0.03);h.param.c4.fill(0.002);
  h.Deltat=1e-7;
  S s;s.param.epsilon=1;s.param.Deltat=h.Deltat;initializeMovingFromHydrostatic(s,h);
  h.solveConductionLAPACKE();s.advanceAcceptedStep();
  for(int i=3;i<20;++i)for(int f=0;f<3;++f)
    check(std::abs(s.value(i,f,S::U)/h.U[f][i]-1)<1e-11,"thermal stencil/legacy mismatch");
}
void parameterRoundTrip() {
  auto s=setup();s.param.q=.123;s.param.c1.fill(.314);s.param.mass={1,2,3};
  s.param.binary_formation=BINARY_FORMATION_POWER_LAW;s.param.capture_coefficient=.42;
  s.param.reflecting_boundary=1;s.param.heating_dispersion=HEATING_RELATIVE_DISPERSION;
  const std::string dir="output/check_moving_parameters/";
  std::filesystem::create_directories(dir);save_param_for_Mathematica(s.param,dir);
  MovingThreeFluidParam result{};std::ifstream file(dir+"param.dat",std::ios::binary);
  file.read(reinterpret_cast<char*>(&result),sizeof(result));
  check(file.gcount()==sizeof(result)&&std::memcmp(&result,&s.param,sizeof(result))==0,"parameter round trip");
}
void singleSplitAgreement() {
  std::array<S,2> sims;
  for(int k=0;k<2;++k) {
    ThreeFluidSim h;h.initSolver(40);h.initPlummer(k?0.5:1,1e-10,k?1:1e-10,1,1);
    h.param.md=h.param.ms;h.param.c2[FS]=h.param.c2[FD]=.17;
    sims[k].param.epsilon=1e-9;sims[k].param.Deltat=1e-5;sims[k].param.max_timestep=1e-5;
    sims[k].param.maxSteps=100;initializeMovingFromHydrostatic(sims[k],h);
    auto observe=[](const S&){};sims[k].evolve(observe);
  }
  for(int i=0;i<sims[0].zones();++i) {
    double total[2]{};
    for(int k=0;k<2;++k)for(int f=0;f<3;++f)total[k]+=sims[k].value(i,f,S::RHO);
    check(std::abs(total[0]/total[1]-1)<1e-8,"single/split total-density equivalence");
    check(std::abs(sims[0].value(i,FS,S::U)/sims[1].value(i,FS,S::U)-1)<1e-8,"single/split dispersion equivalence");
  }
}
void checkerboardDamping() {
  // Subtract an unperturbed twin to exclude the physical open-boundary wave.
  // No conduction is available to hide a defective acoustic coupling.
  constexpr int n=64;constexpr double amplitude=1e-6;
  for(double epsilon:{0.03,1e-9})for(int mode:{S::VEL,S::U}) {
    std::vector<double> faces(n+1),initial(12*n,0);
    for(int i=0;i<=n;++i)faces[i]=double(i)/n;
    for(int i=0;i<n;++i)for(int f=0;f<3;++f) {
      initial[S::index(i,f,S::RHO)]=1e-8;
      initial[S::index(i,f,S::U)]=1;
    }
    const double sound=std::sqrt((10./9.)/epsilon);
    S base;base.param.epsilon=epsilon;base.param.Deltat=10/(n*sound);
    base.param.maxSteps=1;base.initialize(faces,initial,true);S perturbed=base;
    for(int i=0;i<n;++i)for(int f=0;f<3;++f)
      perturbed.state[S::index(i,f,mode)]+=amplitude*(i%2?1:-1)*(mode==S::VEL?sound:1);
    auto obs=[](const S&){};base.evolve(obs);perturbed.evolve(obs);
    double projection=0;
    for(int i=n/4;i<3*n/4;++i) {
      double delta=mode==S::VEL?
        (perturbed.value(i,0,S::VEL)-base.value(i,0,S::VEL))/sound:
        (perturbed.value(i,0,S::RHO)*perturbed.value(i,0,S::U)-
          base.value(i,0,S::RHO)*base.value(i,0,S::U))/1e-8;
      projection+=(i%2?1:-1)*delta;
    }
    const double ratio=std::abs(projection)/(n/2*amplitude);
    std::cout<<"checkerboard epsilon="<<epsilon<<" field="<<mode<<" amplification="<<ratio<<'\n';
    check(ratio<0.1,"undamped acoustic checkerboard");
  }
}
void noConductionEvolution() {
  auto s=setup(80);s.param.epsilon=1e-9;s.param.c1.fill(0);s.param.c2.fill(0);s.param.c4.fill(0);
  s.param.maxTime=0.12;s.param.maxSteps=2000;
  auto observer=[](const S& x) {check(x.energy_ledger_error<1e-10,"adiabatic energy ledger");};
  s.evolve(observer);
  check(s.totalTime>=s.param.maxTime,"no-conduction evolution stalled");
}
void rawMassTransport() {
  for(double velocity:{-.2,.2}) {
    auto s=setup(32);
    for(int i=0;i<s.zones();++i)for(int f=0;f<NF;++f)s.state[S::index(i,f,S::VEL)]=velocity;
    s.assembleStep();
    for(int i=1;i+1<s.zones();++i)for(int f=0;f<NF;++f) {
      // With constant velocity the advective LLF flux is exactly upwind.
      // A frozen primitive offset would change this mass transport.
      const double left=s.value(velocity>0?i-1:i,f,S::RHO)*velocity;
      const double right=s.value(velocity>0?i:i+1,f,S::RHO)*velocity;
      const double expected=-s.Deltat*(std::pow(s.faces()[i+1],2)*right-std::pow(s.faces()[i],2)*left);
      check(std::abs(s.rightHandSide()[S::index(i,f,S::RHO)]-expected)<1e-12*std::abs(expected)+1e-24,
            "reference offsets contaminated mass transport");
    }
  }
}
void contractionRegression(int zones,const std::string& directory) {
  struct Initial {
    long long zones=150;
    double rho0=1,xi1=1e-10,xi2=1,zeta1=1,zeta2=.3;
    double ms=1e-6,mb_over_ms=2,md_over_ms=1e-10;
  } initial;
  initial.zones=zones;check(zones>=16,"too few contraction cells");
  check(!std::filesystem::exists(directory+"/snapshots.csv"),"contraction output exists");
  ThreeFluidSim hydro;hydro.initSolver(zones);
  hydro.param.ms=initial.ms;hydro.param.mb=initial.ms*initial.mb_over_ms;
  hydro.param.md=initial.ms*initial.md_over_ms;
  hydro.initPlummer(initial.rho0,initial.xi1,initial.xi2,initial.zeta1,initial.zeta2);
  double total=0;for(int f=0;f<NF;++f)total+=hydro.Menc[f][zones-1];
  hydro.initCoeffs(total/initial.ms,initial.mb_over_ms,initial.md_over_ms);
  hydro.param.c1.fill(0);hydro.param.c4.fill(0);
  S s;s.param.epsilon=std::pow(3*initial.ms*std::log(0.8*total/(hydro.param.ms+hydro.param.md)),2);
  s.param.maxTime=5;s.param.maxSteps=20000;initializeMovingFromHydrostatic(s,hydro);
  MovingComparisonObserver observer(directory,{});
  std::filesystem::create_directories(directory+"/initialization");
  std::filesystem::create_directories(directory+"/observer");
  save_param_for_Mathematica(initial,directory+"/initialization/");
  save_param_for_Mathematica(observer.param,directory+"/observer/");
  save_param_for_Mathematica(s.param,directory+"/");
  try{s.evolve(observer);}catch(...){observer.saveFinal();throw;}
  observer.saveFinal();
  std::cout<<"contraction zones="<<zones<<" time="<<s.totalTime<<" steps="<<s.step<<'\n';
  check(s.totalTime>=s.param.maxTime,"contraction regression did not finish");
  for(int f=0;f<NF;++f) {
    int turns=0;
    for(int i=1;i+1<s.zones();++i)
      turns+=(s.value(i,f,S::VEL)-s.value(i-1,f,S::VEL))*
        (s.value(i+1,f,S::VEL)-s.value(i,f,S::VEL))<0;
    check(turns<=2,"checkerboard returned in contraction regression");
  }
}
void benchmark() {
  for(int n:{50,100,200,400,800}) {
    auto s=setup(n);s.param.maxSteps=50;s.param.max_timestep=s.param.Deltat;
    auto obs=[](const S&){};
    const auto start=std::chrono::steady_clock::now();s.evolve(obs);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"zones="<<n<<" seconds_per_step="<<seconds/s.step
      <<" band_bytes="<<s.bandMatrix().size()*sizeof(double)<<'\n';
  }
}
void exportFixture(const std::string& directory) {
  std::filesystem::create_directories(directory);
  constexpr int n=5;S s;s.param.epsilon=0.1;s.param.q=0.01;s.param.Deltat=0.01;
  s.param.c2.fill(0.07);s.param.c1.fill(0.003);s.param.c4.fill(0.002);
  s.param.mass={1,2,0.1};
  std::vector<double> faces{0,0.4,0.7,1.1,1.6,2.2},initial(12*n,0);
  for(int i=0;i<n;++i)for(int f=0;f<3;++f) {
    initial[S::index(i,f,S::RHO)]=1+0.1*f+0.04*i;
    initial[S::index(i,f,S::U)]=1+0.2*f+0.03*i;
    initial[S::index(i,f,S::VEL)]=0.2+0.03*f+0.01*i;
    if(i==n-1&&f>0) initial[S::index(i,f,S::VEL)]=(f==1?2:-4)*
      std::sqrt((10./9.)*initial[S::index(i,f,S::U)]/s.param.epsilon);
  }
  s.initialize(faces,initial,false);s.assembleStep();
  write_to_file(s.state,directory+"/state.dat");write_to_file(s.faces(),directory+"/faces.dat");
  const Eigen::MatrixXd a=dense(s);
  write_to_file(a,directory+"/matrix.dat");write_to_file(s.rightHandSide(),directory+"/rhs.dat");
}
void auditHydro() {
  std::cout<<std::setprecision(17);
  for(double dt:{1e-3,1e-4,1e-5,1e-6,1e-7,0.0}) {
    ThreeFluidSim s;MovingComparisonInit p;p.sample=2;initializeMovingComparison(s,p);
    s.param.c2.fill(0);s.Deltat=dt; // isolate local c1 exchange from conduction
    auto moments=[&]() {
      std::array<double,3> result{};
      for(int f=0;f<3;++f)for(int i=0;i<s.param.N;++i) {
        const double lo=i?s.R[f][i-1]:0,hi=s.R[f][i];
        const double mass=s.Rho[f][i]*(hi*hi*hi-lo*lo*lo)/3;
        result[0]+=mass;result[1]+=mass*s.U[f][i];
        result[2]+=mass*std::log(s.P[f][i]/std::pow(s.Rho[f][i],5./3.));
      }return result;
    };
    const auto initial=moments();
    const double u0=s.U[FS][0];
    s.solveConductionLAPACKE();const double heat_u=s.U[FS][0];
    const auto heated=moments();s.projectHydrostatic();const auto projected=moments();
    const double project_u=s.U[FS][0];s.realign();const auto aligned=moments();
    std::cout<<"dt="<<dt<<" align_dM="<<aligned[0]-projected[0]
      <<" align_dE="<<aligned[1]-projected[1]<<" align_dEntropy="<<aligned[2]-projected[2]
      <<" core_dU_heat_absolute="<<heat_u-u0;
    if(dt>0)std::cout<<" conduction_dE_dt="<<(heated[1]-initial[1])/dt
      <<" align_dM_dt="<<(aligned[0]-projected[0])/dt
      <<" align_dE_dt="<<(aligned[1]-projected[1])/dt
      <<" align_dEntropy_dt="<<(aligned[2]-projected[2])/dt
      <<" core_dU_heat="<<(heat_u-u0)/dt<<" core_dU_project="<<(project_u-heat_u)/dt
      <<" core_dU_align="<<(s.U[FS][0]-project_u)/dt;
    std::cout<<'\n';
  }
}
int main(int argc,char**argv) {try {
  if(argc==4&&std::string(argv[1])=="--contraction") {contractionRegression(std::stoi(argv[2]),argv[3]);return 0;}
  if(argc==3&&std::string(argv[1])=="--export") {exportFixture(argv[2]);return 0;}
  if(argc==2&&std::string(argv[1])=="--audit-hydro") {auditHydro();return 0;}
  if(argc==2&&std::string(argv[1])=="--benchmark") {benchmark();return 0;}
  bandSolve();massLedger();
  std::cout<<"checking capture ledger"<<std::endl;captureLedger();
  std::cout<<"checking moving Statler"<<std::endl;movingStatlerSmoke();
  equilibriumAndTide();splitSymmetry();failEarly();dilutedReference();tidalEvolution();
  reflectingAndRelativeHeating();heggieSmoke();
  thermalAgreement();parameterRoundTrip();singleSplitAgreement();checkerboardDamping();noConductionEvolution();rawMassTransport();
  std::cout<<"moving checks passed: band/dense, mass/heat export, equilibrium, q, split, fail-early\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
