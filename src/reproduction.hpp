#pragma once
#include "three_fluid.hpp"

// Declarative example settings; physical metadata are saved in ThreeFluidParam.
ThreeFluidParam statlerParameters(bool direct_heating = true);
void initializeCaptureCluster(ThreeFluidSim& sim, const ThreeFluidParam& settings);
