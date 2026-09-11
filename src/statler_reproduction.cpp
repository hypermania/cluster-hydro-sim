#include "statler_reproduction.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace {
constexpr double solar_radius_pc = 696300.0/3.0856775814913673e13;
constexpr double gravitational_constant = 4.498502151575286e-3; // pc^3 / Msun / Myr^2
constexpr double binary_mass_ratio = 2, dm_mass_ratio = 1e-10;
}

ThreeFluidParam statlerParameters() {
  ThreeFluidParam p;
  p.N=500;p.maxSteps=2000000;p.StopDensity=1e8;
  p.binary_formation=BINARY_FORMATION_POWER_LAW;
  return p;
}

StatlerObserverParam statlerObserverParameters(const StatlerInitParam& p) {
  for(double x:{p.stellar_number,p.stellar_mass_msun,p.stellar_radius_rsun,
                p.length_unit_pc,p.reference_trh_myr})
    if(!std::isfinite(x)||x<=0) throw std::invalid_argument("invalid initializer scales");
  StatlerObserverParam result;
  result.reference_coulomb_log=std::log(0.8*p.stellar_number);
  if(result.reference_coulomb_log<=0) throw std::invalid_argument("nonpositive Coulomb log");
  result.time_unit_myr=std::sqrt(p.stellar_number/3)*std::pow(p.length_unit_pc,1.5)/
    (std::sqrt(gravitational_constant*p.stellar_mass_msun)*result.reference_coulomb_log);
  result.time_unit_over_trh=result.time_unit_myr/p.reference_trh_myr;
  result.snapshot_times_trh[5]=13000./p.reference_trh_myr;
  return result;
}

void initializeCaptureCluster(ThreeFluidSim& s,const StatlerInitParam& init) {
  const auto units=statlerObserverParameters(init);
  auto& p=s.param;
  if(p.N<3) throw std::invalid_argument("invalid grid size");
  s.initSolver(p.N);
  const auto& v=init.plummer;
  s.initPlummer(v[0],v[1],v[2],v[3],v[4]);
  double mass=0;for(int f=0;f<NF;++f)mass+=s.Menc[f][p.N-1];
  p.ms=mass/init.stellar_number;p.mb=binary_mass_ratio*p.ms;p.md=dm_mass_ratio*p.ms;
  const auto formation=p.binary_formation, tidal=p.tidal_cutoff;
  s.initCoeffs(init.stellar_number,binary_mass_ratio,dm_mass_ratio);
  p.binary_formation=formation;p.tidal_cutoff=tidal;
  p.c2[FD]=0;
  for(int f=0;f<NF;++f) {
    p.c1[f*NF+FD]=p.c1[FD*NF+f]=0;
    p.c4[f*NF+FD]=p.c4[FD*NF+f]=0;
  }
  // Maxwellian average of Statler et al. (1987), Eq. (2.1), including the
  // unordered-pair factor 1/2; u_s=3 sigma_s^2/2. See the reproduction notes.
  const double coefficient=107*std::tgamma(0.9)/
    (50*std::pow(2.,0.7)*std::pow(3.,0.4)*std::pow(std::numbers::pi,1.5));
  p.capture_coefficient=coefficient*
    std::pow(init.stellar_radius_rsun*solar_radius_pc/init.length_unit_pc,0.9)*
    std::pow(p.ms,-1.9)/units.reference_coulomb_log;
  s.Deltat=p.Deltat;s.totalTime=0;s.step=0;
  s.cumulative_formed_binaries=0;
  s.validateEvolution();
}
