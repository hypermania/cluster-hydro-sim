#include "three_fluid.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

double ThreeFluidSim::centralDensity() const {
  if (options.central_density_measure == 1) return Rho[FS][0]+Rho[FB][0];
  return std::max({Rho[FS][0],Rho[FB][0],Rho[FD][0]});
}

bool ThreeFluidSim::stopCondition() const {
  return centralDensity()>StopDensity || step>=maxSteps || totalTime>=maxTime;
}

void ThreeFluidSim::sanityCheck() const {
  for(int f=0;f<NF;++f) {
    for(const auto* x:{&R[f],&Rho[f],&U[f],&P[f],&Menc[f]})
      if(!x->allFinite()) throw std::runtime_error("nonfinite simulation state");
    if((R[f].array()<=0).any() || (Rho[f].array()<=0).any() ||
       (U[f].array()<=0).any() || (P[f].array()<=0).any())
      throw std::runtime_error("nonpositive simulation state");
    for(int j=1;j<N;++j) if(R[f][j]<=R[f][j-1])
      throw std::runtime_error("unordered shell radii");
  }
}

void ThreeFluidSim::validateEvolution() const {
  auto positive=[](double x){return std::isfinite(x)&&x>0;};
  if(N<3 || maxSteps<0 || !std::isfinite(maxTime) || maxTime<0 ||
     !std::isfinite(totalTime) || totalTime<0 || !positive(Deltat) ||
     !positive(thres) || !positive(StopDensity) || !positive(options.max_timestep) ||
     options.relaxation_passes<1 || options.max_retries<0 ||
     !positive(options.donor_fraction_limit) || options.donor_fraction_limit>=1 ||
     !positive(options.retry_safety) || options.retry_safety>=1 ||
     !positive(options.density_change_tolerance) || !positive(options.error_control_floor) ||
     !positive(options.timestep_growth_min) || options.timestep_growth_min>1 ||
     !positive(options.timestep_growth_max) || options.timestep_growth_max<1 ||
     options.timestep_controller<0 || options.timestep_controller>1 ||
     options.central_density_measure<0 || options.central_density_measure>1)
    throw std::invalid_argument("invalid evolution controls");
  if(binary_formation==BINARY_FORMATION_POWER_LAW &&
     (!positive(ms)||!positive(mb)||!std::isfinite(options.capture_coefficient)||
      options.capture_coefficient<0||!std::isfinite(options.capture_energy_exponent)||
      !std::isfinite(options.captured_specific_energy_fraction)||
      options.captured_specific_energy_fraction<0||options.captured_specific_energy_fraction>1))
    throw std::invalid_argument("invalid power-law formation settings");
  sanityCheck();
}

double ThreeFluidSim::captureNumberRate(int j) const {
  return options.capture_coefficient*std::pow(Rho[FS][j],2)
    /std::pow(U[FS][j],options.capture_energy_exponent);
}

double ThreeFluidSim::prepareFormationSource() {
  double fraction=0;
  for(int j=0;j<N;++j) {
    formationSource[j]=mb*captureNumberRate(j)*Deltat;
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
  for(int j=0;j<N;++j) {
    const double inner=j==0?0:R[FS][j-1];
    const double volume=(std::pow(R[FS][j],3)-std::pow(inner,3))/3;
    const double transfer=formationSource[j];
    const double birth_u=options.captured_specific_energy_fraction*U[FS][j];
    formed+=transfer*volume/mb;
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
    for(int pass=0;pass<options.relaxation_passes;++pass) solveRelaxationLAPACKE(f);
}

void ThreeFluidSim::selectNextTimestep(double used,double change,double old_density) {
  double factor;
  if(options.timestep_controller==1) {
    const double density_change=std::abs(std::log(centralDensity()/old_density));
    const double control=std::max({change/thres,
      density_change/options.density_change_tolerance,options.error_control_floor});
    factor=std::clamp(1./control,options.timestep_growth_min,options.timestep_growth_max);
    Deltat=std::min(options.max_timestep,used*factor);
  } else {
    // Preserve legacy arithmetic for positive change; zero change has bounded growth.
    Deltat=change>0?used*thres/change:used*options.timestep_growth_max;
    Deltat=std::min(options.max_timestep,Deltat);
  }
  if(!std::isfinite(Deltat)||Deltat<=0) throw std::runtime_error("invalid next timestep");
}

void ThreeFluidSim::advanceAcceptedStep() {
  const double old_density=centralDensity();
  for(int f=0;f<NF;++f){trialU[f]=U[f];trialP[f]=P[f];}
  Deltat=std::min({Deltat,options.max_timestep,maxTime-totalTime});
  if(tidal_cutoff!=TIDAL_CUTOFF_OFF && tidal_cutoff_factor*Deltat>=1)
    throw std::runtime_error("tidal sink would exhaust density");
  double change=0;
  for(long long retry=0;;++retry) {
    if(!std::isfinite(Deltat)||Deltat<=0 || totalTime+Deltat==totalTime)
      throw std::runtime_error("timestep cannot advance time");
    solveConductionLAPACKE();
    sanityCheck();
    change=conductionChange();
    const double fraction=binary_formation==BINARY_FORMATION_POWER_LAW?
      prepareFormationSource():0;
    if(fraction<=options.donor_fraction_limit) break;
    for(int f=0;f<NF;++f){U[f]=trialU[f];P[f]=trialP[f];}
    if(retry>=options.max_retries) throw std::runtime_error("formation retry limit");
    Deltat*=options.retry_safety*options.donor_fraction_limit/fraction;
    ++rejected_steps;
  }
  if(tidal_cutoff!=TIDAL_CUTOFF_OFF) applyTidalCutoff();
  if(binary_formation!=BINARY_FORMATION_OFF) applyBinaryFormation();
  projectHydrostatic();
  realign();
  sanityCheck();
  const double used=Deltat;
  totalTime+=used;
  ++step;
  selectNextTimestep(used,change,old_density);
}
