#pragma once
#include "moving_three_fluid.hpp"
#include <filesystem>

struct HeggieInitParam {
  long long model=0; // 0 Fig.1/2, 1 Fig.3, 2 Fig.4 BS, 3 Fig.4 BS+BB
  long long zones=400;
  double f_star=0.005;
  double binary_number_fraction=0.03;
  double coulomb_log=10;
  double conductivity=0.104;
  double stellar_number=3e5;
  double first_face=1e-6, last_face=1e4;
};
struct HeggieObserverParam {
  double initial_mass=1./3.;
  double time_unit_over_trh=1;
  double length_unit=3*std::numbers::pi/16; // G=M=1, initial E=-1/4 Plummer units
  double density_unit=1;
};

inline double movingHalfMassRadius(const MovingThreeFluidSim& s) {
  double total=0,previous=0;
  for(int f=0;f<NF;++f)total+=s.value(s.zones()-1,f,s.MASS);
  for(int i=0;i<s.zones();++i) {
    double enclosed=0;for(int f=0;f<NF;++f)enclosed+=s.value(i,f,s.MASS);
    if(enclosed>=total/2) {
      const double fraction=(total/2-previous)/(enclosed-previous);
      return std::cbrt(std::pow(s.faces()[i],3)+fraction*
        (std::pow(s.faces()[i+1],3)-std::pow(s.faces()[i],3)));
    }
    previous=enclosed;
  }
  throw std::runtime_error("half-mass radius not found");
}

inline HeggieObserverParam initializeHeggie(MovingThreeFluidSim& s,const HeggieInitParam& p) {
  if(p.model<0||p.model>3||p.zones<3||p.first_face<=0||p.last_face<=p.first_face||
     p.f_star<0||p.coulomb_log<=0||p.conductivity<=0||p.stellar_number<=0||
     p.binary_number_fraction<=0||p.binary_number_fraction>=1)
    throw std::invalid_argument("invalid Heggie initializer");
  for(double x:{p.first_face,p.last_face,p.f_star,p.coulomb_log,p.conductivity,
                p.stellar_number,p.binary_number_fraction})
    if(!std::isfinite(x))throw std::invalid_argument("nonfinite Heggie initializer");
  const double ratio=p.model==0?1e-12:p.model==1?0.01:
                     2*p.binary_number_fraction/(1-p.binary_number_fraction);
  const std::array<double,3> central{1/(1+ratio),ratio/(1+ratio),1e-20};
  std::vector<double> faces(p.zones+1,0),state(12*p.zones,0);
  for(int i=0;i<p.zones;++i) {
    faces[i+1]=p.first_face*std::pow(p.last_face/p.first_face,double(i)/(p.zones-1));
    const double r=(faces[i]+faces[i+1])/2;
    for(int f=0;f<NF;++f) {
      state[s.index(i,f,s.RHO)]=central[f]/std::pow(1+r*r,2.5);
      state[s.index(i,f,s.U)]=1/(12*std::sqrt(1+r*r));
    }
  }
  s.param.reflecting_boundary=1;
  s.param.binary_formation=BINARY_FORMATION_OFF;s.param.capture_coefficient=0;
  s.param.heating_dispersion=HEATING_RELATIVE_DISPERSION;
  s.param.c1.fill(0);s.param.c2.fill(0);s.param.c4.fill(0);
  s.initialize(faces,state,true);
  HeggieObserverParam units;units.initial_mass=0;
  for(int f=0;f<NF;++f)units.initial_mass+=s.value(s.zones()-1,f,s.MASS);
  s.param.mass[FS]=units.initial_mass/p.stellar_number;
  s.param.mass[FB]=(p.model==0?1:2)*s.param.mass[FS];
  s.param.mass[FD]=s.param.mass[FS]; // trace clone of singles, not dark matter
  s.param.epsilon=std::pow(3*s.param.mass[FS]*p.coulomb_log,2);
  for(int f=0;f<NF;++f)s.param.c2[f]=2*std::sqrt(2./3.)*p.conductivity*
    s.param.mass[f]/s.param.mass[FS];
  if(p.model==0) {
    // One-component gas: fixed f* heating, no segregation or formation.
    // The negligible other fluids follow the identical single-fluid solution.
    s.param.c4.fill(3.807*p.f_star/(8*std::numbers::pi));
  } else {
    for(int f=0;f<NF;++f)for(int h=0;h<NF;++h)
      if(f!=h)s.param.c1[3*f+h]=1/std::sqrt(3*std::numbers::pi);
    if(p.model>=2) {
      const double bs=3.807/(16*std::numbers::pi*p.coulomb_log);
      s.param.c4[FS*3+FB]=s.param.c4[FD*3+FB]=2*bs/3;
      s.param.c4[FB*3+FS]=s.param.c4[FB*3+FD]=bs/3;
      if(p.model==3)s.param.c4[FB*3+FB]=5*std::sqrt(3.)/(24*std::numbers::pi*p.coulomb_log);
    }
  }
  const double rh=movingHalfMassRadius(s);
  // Eq.15 uses log10 Lambda, whereas the evolution coefficients use ln Lambda.
  units.time_unit_over_trh=1/(.18*std::log(10.)*std::pow(rh,1.5)*std::sqrt(units.initial_mass));
  units.density_unit=1/(4*std::numbers::pi*std::pow(units.length_unit,3)*units.initial_mass);
  s.validate();return units;
}

struct HeggieObserver {
  HeggieObserverParam param;
  std::ofstream history;
  explicit HeggieObserver(const HeggieObserverParam& p,const std::string& dir)
    :param(p),history(dir+"history.dat",std::ios::binary) {
    history.exceptions(std::ios::badbit|std::ios::failbit);
  }
  void operator()(const MovingThreeFluidSim& s) {
    double mass=0,phi=0,rho=0,thermal=0;
    for(int f=0;f<NF;++f) {
      mass+=s.value(s.zones()-1,f,s.MASS);
      rho+=s.value(0,f,s.RHO);thermal+=s.value(0,f,s.RHO)*s.value(0,f,s.U);
      for(int i=0;i<s.zones();++i)
        phi+=.5*s.value(i,f,s.RHO)*(std::pow(s.faces()[i+1],2)-std::pow(s.faces()[i],2));
    }
    const double rc=std::sqrt(6*thermal/(rho*rho));
    auto relaxation=[&](int f) {return .34*std::pow(2./3.,1.5)*12*std::numbers::pi*
      s.param.mass[FS]/s.param.mass[f]*std::pow(s.value(0,f,s.U),1.5)/s.value(0,f,s.RHO);};
    const std::array<double,13> row{s.totalTime,s.totalTime*param.time_unit_over_trh,
      s.value(0,FS,s.RHO)*param.density_unit,s.value(0,FB,s.RHO)*param.density_unit,
      phi,rc*param.length_unit,movingHalfMassRadius(s)*param.length_unit,
      .517*rho*std::pow(rc,3)/(3*param.initial_mass),mass,
      relaxation(FS),relaxation(FB),s.energy_ledger_error,s.Deltat};
    history.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    if(s.step%10000==0) {
      history.flush();std::cout<<"Heggie step="<<s.step<<" t_trh="<<row[1]
        <<" rho_s="<<row[2]<<" Mc/M="<<row[7]<<std::endl;
    }
  }
  void finish(const std::string& dir,const MovingThreeFluidSim& s,const std::string& status) {
    history.flush();std::ofstream out(dir+"status.txt");
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<status<<"\nsteps="<<s.step<<"\ntime_trh="<<std::setprecision(17)
       <<s.totalTime*param.time_unit_over_trh<<'\n';
  }
};

inline void runMovingHeggie(const std::string& directory,const HeggieInitParam& init) {
  if(directory.empty()||std::filesystem::exists(directory+"/param.dat"))
    throw std::invalid_argument("choose a fresh Heggie output directory");
  MovingThreeFluidSim s;s.param.Deltat=1e-6;s.param.max_timestep=1e4;
  s.param.maxSteps=1000000;s.param.StopDensity=1e9;
  const auto units=initializeHeggie(s,init);
  s.param.maxTime=(init.model==0?1000:150)/units.time_unit_over_trh;
  const std::string dir=directory.back()=='/'?directory:directory+"/";
  std::filesystem::create_directories(dir+"initialization");
  std::filesystem::create_directories(dir+"observer");
  save_param_for_Mathematica(s.param,dir);
  save_param_for_Mathematica(init,dir+"initialization/");
  save_param_for_Mathematica(units,dir+"observer/");
  HeggieObserver observer(units,dir);
  try{s.evolve(observer);}catch(...) {observer.finish(dir,s,"FAILED");throw;}
  const bool density=s.centralDensity()>=s.param.StopDensity;
  observer.finish(dir,s,density?"DENSITY_LIMIT":s.totalTime>=s.param.maxTime?"TIME_LIMIT":"STEP_LIMIT");
  if(!density&&s.totalTime<s.param.maxTime)throw std::runtime_error("Heggie step limit reached");
}
