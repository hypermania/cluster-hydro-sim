#include "../src/three_fluid.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool value,const char* message) {
  if(!value)throw std::runtime_error(message);
}
void close(double a,double b,const char* message) {
  require(std::isfinite(a)&&std::abs(a-b)<2e-13*std::max(1.,std::abs(b)),message);
}
double mass(const ThreeFluidSim& s) {
  double total=0;
  for(int f=0;f<NF;++f)total+=s.Menc[f][s.param.N-1];
  return total;
}
ThreeFluidSim fixture() {
  ThreeFluidSim s;s.initSolver(4);
  for(int f=0;f<NF;++f) {
    s.R[f]<<1.,2.,3.,4.;s.Rho[f].setConstant(f+1.);
    s.U[f].setConstant(f+7.);
    s.P[f]=(2./3.)*s.U[f].array()*s.Rho[f].array();
  }
  s.updateEnclosedMass();s.Deltat=.1;
  s.param.tidal_cutoff=TIDAL_CUTOFF_MASS_SCALED;
  s.param.tidal_radius=2.5;s.param.tidal_cutoff_factor=1.;
  return s;
}
// Independent implementation of the existing fixed-radius drain formula.
void referenceDrain(ThreeFluidSim& s,double radius) {
  for(int f=0;f<NF;++f) {
    for(int i=1;i<s.param.N;++i) {
      if(radius<s.R[f][i-1])s.Rho[f][i]*=1-s.param.tidal_cutoff_factor*s.Deltat;
      else if(radius<s.R[f][i]) {
        const double t=(radius-s.R[f][i-1])/(s.R[f][i]-s.R[f][i-1]);
        s.Rho[f][i]*=1-(1-t)*s.param.tidal_cutoff_factor*s.Deltat;
      }
    }
    s.P[f]=(2./3.)*s.U[f].array()*s.Rho[f].array();
  }
  s.updateEnclosedMass();
}
void sameState(const ThreeFluidSim& a,const ThreeFluidSim& b) {
  for(int f=0;f<NF;++f) for(int i=0;i<a.param.N;++i) {
    close(a.Rho[f][i],b.Rho[f][i],"density mismatch");
    close(a.Menc[f][i],b.Menc[f][i],"enclosed mass mismatch");
    close(a.P[f][i],b.P[f][i],"pressure mismatch");
    require(a.U[f][i]==b.U[f][i]&&a.R[f][i]==b.R[f][i],"sink changed U or R");
  }
}
void operatorChecks() {
  auto s=fixture();const double initial=mass(s);
  auto fixed=s, expected=s;fixed.param.tidal_cutoff=TIDAL_CUTOFF_ON;
  referenceDrain(expected,2.5);fixed.applyTidalCutoff();s.applyTidalCutoff();
  sameState(fixed,expected);sameState(s,expected);
  for(int f=0;f<NF;++f)
    require((fixed.Rho[f].array()==expected.Rho[f].array()).all(),"fixed mode changed arithmetic");
  close(s.tidal_initial_mass,initial,"initial total mass not captured");

  // M/M0=1/8 must halve r_J, including the mass of every component.
  for(int f=0;f<NF;++f)s.Rho[f].setConstant((f+1.)/8.);
  s.updateEnclosedMass();s.totalTime=.1;
  expected=s;referenceDrain(expected,1.25);s.applyTidalCutoff();sameState(s,expected);
  close(s.param.tidal_radius,2.5,"configured initial radius was overwritten");
  close(s.tidal_initial_mass,initial,"reference mass was overwritten");
  // The next call must use the newly depleted mass, without compounding r_J.
  expected=s;referenceDrain(expected,2.5*std::cbrt(mass(s)/initial));
  s.applyTidalCutoff();sameState(s,expected);

  // Mass increase as well as loss follows the cube-root law. Components can
  // change unequally; a single radius must be computed before any is drained.
  s.Rho[FD]*=100.;s.updateEnclosedMass();
  expected=s;referenceDrain(expected,2.5*std::cbrt(mass(s)/initial));
  s.applyTidalCutoff();sameState(s,expected);

  auto off=fixture();off.param.tidal_cutoff=TIDAL_CUTOFF_OFF;
  expected=off;off.applyTidalCutoff();sameState(off,expected);
  require(off.tidal_initial_mass==0,"off mode captured reference mass");
  s.initSolver(4);require(s.tidal_initial_mass==0,"new solver retained old mass");
}
void invalidChecks() {
  for(int kind=0;kind<8;++kind) {
    auto s=fixture();
    switch(kind) {
    case 0:s.param.tidal_cutoff=99;break;
    case 1:s.param.tidal_radius=0;break;
    case 2:s.param.tidal_cutoff_factor=-1;break;
    case 3:s.Deltat=1;break;
    case 4:s.Menc[FD][3]=std::numeric_limits<double>::quiet_NaN();break;
    case 5:s.tidal_initial_mass=-1;break;
    case 6:s.totalTime=1;break; // Lost reference mass on restart.
    case 7:s.Deltat=std::numeric_limits<double>::quiet_NaN();break;
    }
    const auto before=s.Rho;bool failed=false;
    try{s.applyTidalCutoff();}catch(const std::exception&){failed=true;}
    require(failed,"invalid tidal input accepted");
    for(int f=0;f<NF;++f)require((before[f].array()==s.Rho[f].array()).all(),"invalid sink modified density");
  }
}
void evolutionChecks() {
  ThreeFluidSim s;s.initSolver(80);s.initPlummer(.5,1e-10,1.,1.,1.);
  s.initCoeffs(1e6,2.,.1);s.param.c4.fill(0);
  s.param.md=.1*s.param.ms;
  s.param.tidal_cutoff=TIDAL_CUTOFF_MASS_SCALED;
  s.param.tidal_radius=2.;s.param.tidal_cutoff_factor=10;
  s.param.maxSteps=30;s.param.max_timestep=.001;
  const double initial=mass(s);int calls=0;
  auto observer=[&](const ThreeFluidSim& state){
    ++calls;
    if(state.step>0)close(state.tidal_initial_mass,initial,"evolve changed reference mass");
  };
  s.evolve(observer);s.sanityCheck();
  require(calls==31&&s.step==30,"unexpected evolve history");
  require(mass(s)<initial,"tidal evolve did not lose mass");
  // Continuing the same object must retain the original normalization.
  s.evolve(observer);close(s.tidal_initial_mass,initial,"continuation reset reference mass");
  std::filesystem::create_directories("output/check_tidal_cutoff_data");
  s.saveParams("output/check_tidal_cutoff_data/");
  const auto restored=read_from_file<std::vector<ThreeFluidParam>>("output/check_tidal_cutoff_data/param.dat");
  require(restored&&restored->size()==1,"parameter round trip failed");
  require((*restored)[0].tidal_cutoff==TIDAL_CUTOFF_MASS_SCALED,"mode not serialized");
  close((*restored)[0].tidal_radius,2.,"initial radius not serialized");
}
int main() {
  try{operatorChecks();invalidChecks();evolutionChecks();
    std::cout<<"tidal cutoff checks passed: fixed, cube-root, repeated steps, total mass, validation, evolve, serialization\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
