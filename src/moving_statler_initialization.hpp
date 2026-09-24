#pragma once
#include "moving_initialization.hpp"
#include "statler_reproduction.hpp"

// Statler's radial domain (Section III d), independent of evolution controls.
struct MovingStatlerGrid {
  long long zones=500;
  double first_face=1e-6, last_face=1e4;
};

inline void initializeMovingCaptureCluster(MovingThreeFluidSim& sim,
    const StatlerInitParam& initial,const MovingStatlerGrid& grid) {
  if(grid.zones<3||!std::isfinite(grid.first_face)||!std::isfinite(grid.last_face)||
     grid.first_face<=0||grid.last_face<=grid.first_face)
    throw std::invalid_argument("invalid moving Statler grid");
  if(initial.plummer[3]!=1||initial.plummer[4]!=1)
    throw std::invalid_argument("analytic moving Plummer initializer requires common scale radii");
  ThreeFluidSim hydro;hydro.param=statlerParameters();hydro.param.N=grid.zones;
  initializeCaptureCluster(hydro,initial);
  initializeMovingFromHydrostatic(sim,hydro);
  const auto observing=statlerObserverParameters(initial);
  sim.param.epsilon=std::pow(3*hydro.param.ms*observing.reference_coulomb_log,2);
  std::vector<double> faces(grid.zones+1,0),state(12*grid.zones,0);
  const double total_density=initial.plummer[0]*(1+initial.plummer[1]+initial.plummer[2]);
  for(int i=0;i<grid.zones;++i) {
    faces[i+1]=grid.first_face*std::pow(grid.last_face/grid.first_face,double(i)/(grid.zones-1));
    const double r=(faces[i]+faces[i+1])/2;
    for(int f=0;f<NF;++f) {
      const double fraction=f==FS?1:initial.plummer[f];
      state[sim.index(i,f,sim.RHO)]=initial.plummer[0]*fraction/std::pow(1+r*r,2.5);
      // Analytic isotropic Plummer pressure, U=3 sigma_1D^2/2.
      state[sim.index(i,f,sim.U)]=total_density/(12*std::sqrt(1+r*r));
    }
  }
  // Preserve the prescribed stellar number on this finite-volume grid. The
  // formation coefficient scales as m_s^-1.9 in the same external conversion.
  double total_mass=0;
  for(int i=0;i<grid.zones;++i)for(int f=0;f<NF;++f)
    total_mass+=(std::pow(faces[i+1],3)-std::pow(faces[i],3))/3*state[sim.index(i,f,sim.RHO)];
  const double ratio=(total_mass/initial.stellar_number)/sim.param.mass[FS];
  for(double& m:sim.param.mass)m*=ratio;
  sim.param.capture_coefficient*=std::pow(ratio,-1.9);
  sim.param.epsilon*=ratio*ratio;
  sim.initialize(faces,state,true);
}
