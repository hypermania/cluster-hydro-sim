#include "three_fluid.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

double ThreeFluidSim::centralDensity() const {
  return std::max({Rho[FS][0],Rho[FB][0],Rho[FD][0]});
}

bool ThreeFluidSim::stopCondition() const {
  return centralDensity()>param.StopDensity || step>=param.maxSteps || totalTime>=param.maxTime;
}

void ThreeFluidSim::sanityCheck() const {
  for(int f=0;f<NF;++f) {
    for(const auto* x:{&R[f],&Rho[f],&U[f],&P[f],&Menc[f]})
      if(!x->allFinite()) throw std::runtime_error("nonfinite simulation state");
    if((R[f].array()<=0).any() || (Rho[f].array()<=0).any() ||
       (U[f].array()<=0).any() || (P[f].array()<=0).any())
      throw std::runtime_error("nonpositive simulation state");
    for(int j=1;j<param.N;++j) if(R[f][j]<=R[f][j-1])
      throw std::runtime_error("unordered shell radii");
  }
}

void ThreeFluidSim::validateEvolution() const {
  // Always enforce storage compatibility, even when full-state scans are off.
  for(int f=0;f<NF;++f)
    for(const auto* x:{&R[f],&Rho[f],&U[f],&P[f],&Menc[f]})
      if(x->size()!=param.N) throw std::invalid_argument("initialize solver after changing N");
  auto positive=[](double x){return std::isfinite(x)&&x>0;};
  if(param.N<3 || param.maxSteps<0 || !std::isfinite(param.maxTime) || param.maxTime<0 ||
     !std::isfinite(totalTime) || totalTime<0 || !positive(Deltat) ||
     !positive(param.u_change_tolerance) || !positive(param.StopDensity) || !positive(param.max_timestep) ||
     (param.runtime_validation!=0 && param.runtime_validation!=1))
    throw std::invalid_argument("invalid evolution controls");
  if(param.binary_formation==BINARY_FORMATION_POWER_LAW &&
     (!positive(param.ms)||!positive(param.mb)||!std::isfinite(param.capture_coefficient)||
      param.capture_coefficient<0))
    throw std::invalid_argument("invalid power-law formation settings");
  if(param.runtime_validation) sanityCheck();
}

double ThreeFluidSim::captureNumberRate(int j) const {
  return param.capture_coefficient*std::pow(Rho[FS][j],2)
    /std::pow(U[FS][j],0.6);
}

void ThreeFluidSim::applyPowerLawFormation() {
  // Compute once per zone, after any density sink. On failure terminate the
  // run: earlier zones may already be modified and must not be saved/resumed.
  double formed=0;
  for(int j=0;j<param.N;++j) {
    const double inner=j==0?0:R[FS][j-1];
    const double volume=(std::pow(R[FS][j],3)-std::pow(inner,3))/3;
    const double transfer=param.mb*captureNumberRate(j)*Deltat;
    if(!std::isfinite(transfer) || transfer<0 || !(transfer<Rho[FS][j]))
      throw std::runtime_error("invalid capture transfer in zone "+std::to_string(j));
    const double birth_u=0.5*U[FS][j];
    formed+=transfer*volume/param.mb;
    U[FB][j]=(Rho[FB][j]*U[FB][j]+transfer*birth_u)/(Rho[FB][j]+transfer);
    Rho[FS][j]-=transfer;Rho[FB][j]+=transfer;
  }
  cumulative_formed_binaries+=formed;
  updateEnclosedMass();
  for(int f=0;f<NF;++f) P[f].array()=(2./3.)*Rho[f].array()*U[f].array();
}

double ThreeFluidSim::conductionChange() const {
  double result=0;
  for(int f=0;f<NF;++f) for(int j=0;j<param.N;++j) {
    const double change=std::abs(U[f][j]-previousU[f][j])/previousU[f][j];
    if(!std::isfinite(change) || previousU[f][j]<=0 || U[f][j]<=0)
      throw std::runtime_error("invalid conduction energy change");
    result=std::max(result,change);
  }
  return result;
}

void ThreeFluidSim::projectHydrostatic() {
  for(int f=0;f<NF;++f)
    for(int pass=0;pass<2;++pass) solveRelaxationLAPACKE(f);
}

void ThreeFluidSim::selectNextTimestep(double used,double change) {
  // if(!std::isfinite(change)||change<0)
  //   throw std::runtime_error("invalid timestep error estimate");
  Deltat=change>0?std::min(param.max_timestep,used*param.u_change_tolerance/change)
                :param.max_timestep;
  if(!std::isfinite(Deltat)||Deltat<=0) throw std::runtime_error("invalid next timestep");
}

void ThreeFluidSim::advanceAcceptedStep() {
  if(!std::isfinite(Deltat)||Deltat<=0)
    throw std::runtime_error("invalid timestep");
  // Deltat=std::min({Deltat,param.max_timestep,param.maxTime-totalTime});
  // if(Deltat<=0 || totalTime+Deltat==totalTime)
  //   throw std::runtime_error("timestep cannot advance time");
  if(param.tidal_cutoff!=TIDAL_CUTOFF_OFF && param.tidal_cutoff_factor*Deltat>=1)
    throw std::runtime_error("tidal sink would exhaust density");
  previousU=U;
  solveConductionLAPACKE();
  const double change=conductionChange();
  if(param.tidal_cutoff!=TIDAL_CUTOFF_OFF) applyTidalCutoff();
  if(param.binary_formation!=BINARY_FORMATION_OFF) applyBinaryFormation();
  projectHydrostatic();
  realign();
  if(param.runtime_validation) sanityCheck();
  const double used=Deltat;
  totalTime+=used;
  ++step;
  selectNextTimestep(used,change);
}
