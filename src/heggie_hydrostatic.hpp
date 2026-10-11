#pragma once
#include "three_fluid.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <stdexcept>

// Initializer/observer settings are separate from dimensionless evolution.
struct HydroHeggieInitParam {
  long long model=0; // 0: Fig.1/2; 1: Fig.3; 2: Fig.4 BS; 3: Fig.4 BS+BB
  long long zones=400;
  double f_star=.005;
  double binary_number_fraction=.03;
  double coulomb_log=10;
  double conductivity=.104;
  double stellar_number=3e5;
};
struct HydroHeggieObserverParam {
  double initial_mass=1./3.;
  double time_unit_over_trh=1;
  double length_unit=3*std::numbers::pi/16;
  double density_unit=1;
  long long primary_component=FS;
};

// Initializer/observer helpers only: ThreeFluidSim evolution is unchanged.
inline double hydroHalfMassRadius(const ThreeFluidSim& s) {
  double total=0,previous=0;
  for(int f=0;f<NF;++f)total+=s.Menc[f][s.param.N-1];
  for(int i=0;i<s.param.N;++i) {
    double enclosed=0;for(int f=0;f<NF;++f)enclosed+=s.Menc[f][i];
    if(enclosed>=total/2) {
      const double inner=i?s.R[FS][i-1]:0;
      return std::cbrt(std::pow(inner,3)+(total/2-previous)/(enclosed-previous)*
                      (std::pow(s.R[FS][i],3)-std::pow(inner,3)));
    }
    previous=enclosed;
  }
  throw std::runtime_error("half-mass radius not found");
}

struct HeggieHydroProfileParam {
  double single_scale=1;
  double binary_scale_ratio=1;
};

inline HydroHeggieObserverParam initializeHydroHeggie(ThreeFluidSim& s,const HydroHeggieInitParam& p,
                                               const HeggieHydroProfileParam& profile={}) {
  const double a=profile.single_scale,b=profile.binary_scale_ratio;
  if(!std::isfinite(a)||!std::isfinite(b)||a<=0||b<=0)
    throw std::invalid_argument("Plummer profile scales must be finite and positive");
  if(p.model<0||p.model>3||p.zones<3||p.zones>100000 ||
     !std::isfinite(p.f_star)||p.f_star<0||
     !std::isfinite(p.binary_number_fraction)||p.binary_number_fraction<=0||p.binary_number_fraction>=1||
     !std::isfinite(p.coulomb_log)||p.coulomb_log<=0||
     !std::isfinite(p.conductivity)||p.conductivity<=0||
     !std::isfinite(p.stellar_number)||p.stellar_number<=0)
    throw std::invalid_argument("invalid Heggie initialization");
  if(p.model<2 && (a!=1 || b!=1))
    throw std::invalid_argument("profile adjustment is supported for Fig.4 only");
  HydroHeggieObserverParam units;
  s.initSolver(p.zones);
  const double ratio=p.model==0?1e-12:p.model==1?.01:2*p.binary_number_fraction/(1-p.binary_number_fraction);
  constexpr double trace=1e-12;
  if(p.model==0) {
    // The stock self-heating kernel is in slot FB. This is ONE effective
    // fluid, not a cluster made entirely of physical binaries.
    s.initPlummer(trace/(1+2*trace),1/trace,1,1,1);
    units.primary_component=FB;
  } else {
    s.initPlummer(1/(1+ratio),ratio/std::pow(b,3),1e-20*(1+ratio),b,1);
  }
  // Fixed-mass hydrostatic homology; retain the native grid and its boundary.
  // The initial grid becomes [0.01*a,1000*a], not the experimental wide grid.
  for(int f=0;f<NF;++f) {
    s.R[f]*=a;s.Rho[f]/=std::pow(a,3);s.P[f]/=std::pow(a,4);s.U[f]/=a;
  }
  s.param.c1.fill(0);s.param.c2.fill(0);s.param.c4.fill(0);
  units.initial_mass=0;
  for(int f=0;f<NF;++f)units.initial_mass+=s.Menc[f][p.zones-1];
  s.param.ms=units.initial_mass/p.stellar_number;
  s.param.mb=(p.model==0?1:2)*s.param.ms;s.param.md=s.param.ms;
  for(int f=0;f<NF;++f)
    s.param.c2[f]=2*std::sqrt(2./3.)*p.conductivity*(f==FB?s.param.mb/s.param.ms:1);
  if(p.model==0) {
    const double single=3.807*p.f_star/(8*std::numbers::pi*std::sqrt(2.))*(1+2*trace);
    for(int f=0;f<NF;++f)s.param.c4[f*NF+FB]=single;
  } else {
    for(int f=0;f<NF;++f)for(int j=0;j<NF;++j)
      if(f!=j)s.param.c1[f*NF+j]=1/std::sqrt(3*std::numbers::pi);
    if(p.model>=2) {
      // Existing kernel: BS denominator sqrt(U_s), BB denominator sqrt(U_b).
      // Match Heggie's relative-dispersion rates at U_s=U_b only.
      const double bs=3.807/(16*std::numbers::pi*p.coulomb_log*std::sqrt(2.));
      s.param.c4[FS*NF+FB]=s.param.c4[FD*NF+FB]=2*bs/3;
      s.param.c4[FB*NF+FS]=s.param.c4[FB*NF+FD]=bs/3;
      if(p.model==3)s.param.c4[FB*NF+FB]=5*std::sqrt(3.)/
        (24*std::numbers::pi*p.coulomb_log*std::sqrt(2.));
    }
  }
  s.param.binary_formation=BINARY_FORMATION_OFF;s.param.tidal_cutoff=TIDAL_CUTOFF_OFF;
  const double rh=hydroHalfMassRadius(s);
  units.time_unit_over_trh=1/(.18*std::log(10.)*std::pow(rh,1.5)*std::sqrt(units.initial_mass));
  units.density_unit=1/(4*std::numbers::pi*std::pow(units.length_unit,3)*units.initial_mass);
  return units;
}

struct HydroHeggieObserver {
  HydroHeggieObserverParam param;
  std::ofstream history, diagnostics;
  explicit HydroHeggieObserver(const HydroHeggieObserverParam& p,const std::string& dir)
    :param(p),history(dir+"history.dat",std::ios::binary),diagnostics(dir+"diagnostics.dat",std::ios::binary) {
    history.exceptions(std::ios::badbit|std::ios::failbit);
    diagnostics.exceptions(std::ios::badbit|std::ios::failbit);
  }
  void operator()(const ThreeFluidSim& s) {
    double mass=0,phi=0,rho=0,thermal=0,energy=0,outer_mass=0;
    for(int f=0;f<NF;++f) {
      mass+=s.Menc[f][s.param.N-1];rho+=s.Rho[f][0];thermal+=s.Rho[f][0]*s.U[f][0];
      for(int i=0;i<s.param.N;++i) {
        const double lo=i?s.R[f][i-1]:0,hi=s.R[f][i];
        phi+=.5*s.Rho[f][i]*(hi*hi-lo*lo);
        energy+=s.Rho[f][i]*s.U[f][i]*(std::pow(hi,3)-std::pow(lo,3))/3;
      }
    }
    // Exact self-energy integral for piecewise constant total shell density.
    double enclosed=0;
    for(int i=0;i<s.param.N;++i) {
      const double lo=i?s.R[FS][i-1]:0,hi=s.R[FS][i];
      double density=0;for(int f=0;f<NF;++f)density+=s.Rho[f][i];
      if(hi>100)outer_mass+=density*(std::pow(hi,3)-std::pow(std::max(lo,100.),3))/3;
      energy-=density*((enclosed-density*std::pow(lo,3)/3)*(hi*hi-lo*lo)/2+
                       density*(std::pow(hi,5)-std::pow(lo,5))/15);
      enclosed+=density*(std::pow(hi,3)-std::pow(lo,3))/3;
    }
    const double rc=std::sqrt(6*thermal/(rho*rho));
    auto relaxation=[&](int f) {return .34*std::pow(2./3.,1.5)*12*std::numbers::pi*
      s.param.ms/(f==FS?s.param.ms:f==FB?s.param.mb:s.param.md)*
      std::pow(s.U[f][0],1.5)/s.Rho[f][0];};
    // Column 11 is total mechanical energy, NOT the moving solver's ledger error.
    const std::array<double,13> row{s.totalTime,s.totalTime*param.time_unit_over_trh,
      s.Rho[param.primary_component][0]*param.density_unit,s.Rho[param.primary_component==FS?FB:FS][0]*param.density_unit,
      phi,rc*param.length_unit,hydroHalfMassRadius(s)*param.length_unit,
      .517*rho*std::pow(rc,3)/(3*param.initial_mass),mass,
      relaxation(FS),relaxation(FB),energy,s.Deltat};
    history.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    const std::array<double,8> extra{row[1],s.U[FS][0]/s.U[FB][0],
      s.R[FS][s.param.N-1],s.R[FB][s.param.N-1],s.R[FD][s.param.N-1],
      s.U[FS][s.param.N-1],s.U[FB][s.param.N-1],outer_mass/mass};
    diagnostics.write(reinterpret_cast<const char*>(extra.data()),sizeof(extra));
    if(s.step%1000==0) {
      history.flush();diagnostics.flush();std::cout<<"Hydro Heggie step="<<s.step<<" t_trh="<<row[1]
        <<" rho_s="<<row[2]<<" mass_ratio="<<mass/param.initial_mass<<" dt="<<s.Deltat<<std::endl;
    }
  }
  void finish(const std::string& dir,const ThreeFluidSim& s,const std::string& status) {
    history.flush();diagnostics.flush();std::ofstream out(dir+"status.txt");
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<status<<"\nsteps="<<s.step<<"\ntime_trh="<<std::setprecision(17)
       <<s.totalTime*param.time_unit_over_trh<<'\n';
  }
};

inline void runHydroHeggie(const std::string& directory,const HydroHeggieInitParam& init,
                          double final_trh=0,double tolerance=1e-3,long long steps=1000000,
                          const HeggieHydroProfileParam& profile={}) {
  if(directory.empty()||std::filesystem::exists(directory))
    throw std::invalid_argument("choose a fresh Heggie output directory");
  if(!std::isfinite(final_trh)||final_trh<0||!std::isfinite(tolerance)||tolerance<=0||steps<1)
    throw std::invalid_argument("invalid Heggie run controls");
  ThreeFluidSim s;
  auto units=initializeHydroHeggie(s,init,profile);
  s.param.Deltat=1e-7;s.param.max_timestep=1e4;s.param.maxSteps=steps;
  s.param.StopDensity=1e9;s.param.u_change_tolerance=tolerance;
  s.param.maxTime=(final_trh>0?final_trh:init.model==0?1000:150)/units.time_unit_over_trh;
  const std::string dir=directory.back()=='/'?directory:directory+"/";
  std::filesystem::create_directories(dir+"initialization");
  std::filesystem::create_directories(dir+"observer");
  std::filesystem::create_directories(dir+"profile");
  save_param_for_Mathematica(s.param,dir);
  save_param_for_Mathematica(init,dir+"initialization/");
  save_param_for_Mathematica(profile,dir+"profile/");
  save_param_for_Mathematica(units,dir+"observer/");
  std::ofstream(dir+"solver.txt")<<"hydrostatic-stock\ncolumn_11=total_mechanical_energy_code\nheating=stock_equal_dispersion_matched\n";
  HydroHeggieObserver observer(units,dir);
  try{s.evolve(observer);}catch(...) {observer.finish(dir,s,"FAILED");throw;}
  const bool density=s.centralDensity()>=s.param.StopDensity;
  observer.finish(dir,s,density?"DENSITY_LIMIT":s.totalTime>=s.param.maxTime?"TIME_LIMIT":"STEP_LIMIT");
  if(!density&&s.totalTime<s.param.maxTime)throw std::runtime_error("Heggie step limit reached");
}
