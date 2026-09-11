#pragma once
#include "three_fluid.hpp"
#include "statler_observer.hpp"

// Dimensional inputs belong to the external initializer, never ThreeFluidSim.
struct StatlerInitParam {
  double stellar_number = 3e5;
  double stellar_mass_msun = 0.7;
  double stellar_radius_rsun = 0.57;
  double length_unit_pc = 1.13;
  double reference_trh_myr = 225;
  std::array<double, 5> plummer = {1, 1e-12, 1e-20, 1, 1};
};

ThreeFluidParam statlerParameters();
StatlerObserverParam statlerObserverParameters(const StatlerInitParam&);
void initializeCaptureCluster(ThreeFluidSim&, const StatlerInitParam&);
