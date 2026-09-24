#pragma once
#include "moving_three_fluid.hpp"

// Orthogonal adapter: copy coefficients and profiles, not run controls or IO.
// epsilon remains explicitly chosen in the external initializer/runner.
inline void initializeMovingFromHydrostatic(MovingThreeFluidSim& moving,
                                           const ThreeFluidSim& hydro) {
  if((hydro.param.binary_formation!=BINARY_FORMATION_OFF&&
      hydro.param.binary_formation!=BINARY_FORMATION_POWER_LAW)||hydro.param.tidal_cutoff!=TIDAL_CUTOFF_OFF)
    throw std::invalid_argument("moving solver supports only OFF/POWER_LAW formation and no stripping");
  moving.param.mass={hydro.param.ms,hydro.param.mb,hydro.param.md};
  moving.param.c1=hydro.param.c1;moving.param.c2=hydro.param.c2;moving.param.c4=hydro.param.c4;
  moving.param.binary_formation=hydro.param.binary_formation;
  moving.param.capture_coefficient=hydro.param.capture_coefficient;
  const int n=hydro.param.N;
  std::vector<double> faces(n+1,0),state(12*n,0);
  for(int i=0;i<n;++i) {
    faces[i+1]=hydro.R[FS][i];
    for(int f=0;f<NF;++f) {
      if(hydro.R[f][i]!=faces[i+1])throw std::invalid_argument("moving solver requires common initial grid");
      state[MovingThreeFluidSim::index(i,f,MovingThreeFluidSim::RHO)]=hydro.Rho[f][i];
      state[MovingThreeFluidSim::index(i,f,MovingThreeFluidSim::U)]=hydro.U[f][i];
    }
  }
  moving.initialize(faces,state,true);
}
