#include "three_fluid.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
constexpr double retry_safety = 0.8;
constexpr int max_retries = 100;
constexpr double growth_min = 0.5, growth_max = 1.5, error_floor = 1e-12;
}

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
     !positive(param.thres) || !positive(param.StopDensity) || !positive(param.max_timestep) ||
     !positive(param.donor_fraction_limit) || param.donor_fraction_limit>=1 ||
     !positive(param.density_change_tolerance) ||
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

double ThreeFluidSim::prepareFormationSource() {
  double fraction=0;
  for(int j=0;j<param.N;++j) {
    formationSource[j]=param.mb*captureNumberRate(j)*Deltat;
    if(!std::isfinite(formationSource[j]) || formationSource[j]<0)
      throw std::runtime_error("invalid formation source");
    fraction=std::max(fraction,formationSource[j]/Rho[FS][j]);
  }
  return fraction;
}

void ThreeFluidSim::applyPowerLawFormation() {
  // Re-evaluate after any preceding density sink. No changes before validation.
  if(prepareFormationSource()>=1)
    throw std::runtime_error("formation would exhaust donor density");
  double formed=0;
  for(int j=0;j<param.N;++j) {
    const double inner=j==0?0:R[FS][j-1];
    const double volume=(std::pow(R[FS][j],3)-std::pow(inner,3))/3;
    const double transfer=formationSource[j];
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
  for(int f=0;f<NF;++f)
    result=std::max(result,((U[f].array()-trialU[f].array()).abs()/trialU[f].array()).maxCoeff());
  if(!std::isfinite(result)) throw std::runtime_error("nonfinite conduction change");
  return result;
}

void ThreeFluidSim::projectHydrostatic() {
  for(int f=0;f<NF;++f)
    for(int pass=0;pass<2;++pass) solveRelaxationLAPACKE(f);
}

void ThreeFluidSim::selectNextTimestep(double used,double change,double old_density) {
  const double density_change=std::abs(std::log(centralDensity()/old_density));
  const double control=std::max({change/param.thres,
    density_change/param.density_change_tolerance,error_floor});
  Deltat=std::min(param.max_timestep,
                 used*std::clamp(1./control,growth_min,growth_max));
  if(!std::isfinite(Deltat)||Deltat<=0) throw std::runtime_error("invalid next timestep");
}

void ThreeFluidSim::advanceAcceptedStep() {
  const double old_density=centralDensity();
  for(int f=0;f<NF;++f){trialU[f]=U[f];trialP[f]=P[f];}
  Deltat=std::min({Deltat,param.max_timestep,param.maxTime-totalTime});
  if(param.tidal_cutoff!=TIDAL_CUTOFF_OFF && param.tidal_cutoff_factor*Deltat>=1)
    throw std::runtime_error("tidal sink would exhaust density");
  double change=0;
  for(long long retry=0;;++retry) {
    if(!std::isfinite(Deltat)||Deltat<=0 || totalTime+Deltat==totalTime)
      throw std::runtime_error("timestep cannot advance time");
    solveConductionLAPACKE();
    if(param.runtime_validation) sanityCheck();
    change=conductionChange();
    const double fraction=param.binary_formation==BINARY_FORMATION_POWER_LAW?
      prepareFormationSource():0;
    if(fraction<=param.donor_fraction_limit) break;
    for(int f=0;f<NF;++f){U[f]=trialU[f];P[f]=trialP[f];}
    if(retry>=max_retries) throw std::runtime_error("formation retry limit");
    Deltat*=retry_safety*param.donor_fraction_limit/fraction;
    ++rejected_steps;
  }
  if(param.tidal_cutoff!=TIDAL_CUTOFF_OFF) applyTidalCutoff();
  if(param.binary_formation!=BINARY_FORMATION_OFF) applyBinaryFormation();
  projectHydrostatic();
  realign();
  if(param.runtime_validation) sanityCheck();
  const double used=Deltat;
  totalTime+=used;
  ++step;
  selectNextTimestep(used,change,old_density);
}
