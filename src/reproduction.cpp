#include "reproduction.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

ThreeFluidParam statlerParameters(bool direct) {
  ThreeFluidParam p{};
  p.N=500;p.initial_profile_values={1,1e-12,1e-20,1,1};
  p.stellar_number=3e5;p.stellar_mass_msun=0.7;p.stellar_radius_rsun=0.57;
  p.length_unit_pc=1.13;p.reference_trh_myr=225;
  p.solar_radius_pc=696300.0/3.0856775814913673e13;
  p.gravitational_constant_pc3_msun_myr2=4.498502151575286e-3;
  p.reference_coulomb_log=std::log(0.8*p.stellar_number);
  p.final_time_trh=1e4;p.maxSteps=2000000;p.StopDensity=1e8;p.thres=1e-3;
  p.Deltat=p.initial_timestep=1e-3;
  p.binary_formation=BINARY_FORMATION_POWER_LAW;p.tidal_cutoff=TIDAL_CUTOFF_OFF;
  p.timestep_controller=1;p.central_density_measure=1;p.direct_heating=direct;
  p.statler_observer=1;p.snapshot_count=9;
  p.snapshot_times_trh={0,1,5,10,20,13000./p.reference_trh_myr,100,1000,10000};
  return p;
}

void initializeCaptureCluster(ThreeFluidSim& s,const ThreeFluidParam& p) {
  for(double x:{p.stellar_number,p.stellar_mass_msun,p.stellar_radius_rsun,
      p.length_unit_pc,p.reference_trh_myr,p.solar_radius_pc,
      p.gravitational_constant_pc3_msun_myr2,p.reference_coulomb_log})
    if(!std::isfinite(x)||x<=0) throw std::invalid_argument("invalid physical reference metadata");
  if(p.N<3 || p.initial_profile!=0 || p.final_time_trh<0 || !std::isfinite(p.final_time_trh))
    throw std::invalid_argument("invalid capture cluster setup");
  s.options=p;
  s.initSolver(p.N);
  const auto& v=p.initial_profile_values;
  s.initPlummer(v[0],v[1],v[2],v[3],v[4]);
  double mass=0;for(int f=0;f<NF;++f)mass+=s.Menc[f][s.N-1];
  s.ms=mass/p.stellar_number;s.mb=2*s.ms;s.md=1e-10*s.ms;
  s.initCoeffs(p.stellar_number,2,1e-10);
  s.binary_formation=p.binary_formation;s.tidal_cutoff=p.tidal_cutoff;
  s.tidal_radius=p.tidal_radius;s.tidal_cutoff_factor=p.tidal_cutoff_factor;
  s.c2[FD]=0;
  for(int f=0;f<NF;++f) {
    s.c1[f][FD]=s.c1[FD][f]=0;s.c4[f][FD]=s.c4[FD][f]=0;
    if(!p.direct_heating)for(int g=0;g<NF;++g)s.c4[f][g]=0;
  }
  const double coefficient=107*std::tgamma(0.9)/
    (50*std::pow(2.,0.7)*std::pow(3.,0.4)*std::pow(std::numbers::pi,1.5));
  s.options.capture_coefficient=coefficient*
    std::pow(p.stellar_radius_rsun*p.solar_radius_pc/p.length_unit_pc,0.9)*
    std::pow(s.ms,-1.9)/p.reference_coulomb_log;
  s.options.time_unit_myr=std::sqrt(p.stellar_number/3)*std::pow(p.length_unit_pc,1.5)/
    (std::sqrt(p.gravitational_constant_pc3_msun_myr2*p.stellar_mass_msun)*p.reference_coulomb_log);
  s.options.time_unit_over_trh=s.options.time_unit_myr/p.reference_trh_myr;
  s.maxTime=p.final_time_trh/s.options.time_unit_over_trh;
  s.maxSteps=p.maxSteps;s.StopDensity=p.StopDensity;s.thres=p.thres;
  s.Deltat=p.initial_timestep;s.totalTime=0;s.step=0;
  s.cumulative_formed_binaries=0;s.rejected_steps=0;
  s.validateEvolution();
}
