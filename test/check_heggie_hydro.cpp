#include "../src/heggie_hydrostatic.hpp"
#include <iostream>

void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}
int main() {try {
  for(int model=0;model<4;++model) {
    HydroHeggieInitParam p;p.zones=160;p.model=model;
    ThreeFluidSim s;auto units=initializeHydroHeggie(s,p);
    require(s.param.binary_formation==BINARY_FORMATION_OFF,"formation must be off");
    require(s.param.tidal_cutoff==TIDAL_CUTOFF_OFF,"tidal removal must be off");
    require(std::abs(units.initial_mass-1./3.)<.002,"initial mass");
    require(std::abs(s.R[FS][0]-.01)<1e-14 && std::abs(s.R[FS][159]-1000)<1e-10,"native grid changed");
    const auto radius=s.R;s.projectHydrostatic();
    for(int f=0;f<NF;++f)require((s.R[f]-radius[f]).cwiseAbs().maxCoeff()<1e-8,"not hydrostatic");
    s.param.maxSteps=20;s.param.maxTime=1;s.param.Deltat=1e-7;
    int count=0;auto obs=[&](const ThreeFluidSim& state){state.sanityCheck();++count;};
    s.evolve(obs);require(count==21,"evolve observer count");
    if(model==0) {
      require(units.primary_component==FB,"effective fluid slot");
      require((s.U[FS].array()/s.U[FB].array()-1).abs().maxCoeff()<1e-8,"trace not a clone");
    }
  }
  HydroHeggieInitParam p;p.zones=200;p.model=3;
  ThreeFluidSim s;const auto unit=initializeHydroHeggie(s,p,{.72,1.5});
  require(std::abs(s.Menc[FB][199]/s.Menc[FS][199]/(.06/.97)-1)<3e-6,"profile changes mass fraction");
  require(std::abs(s.R[FS][0]-.0072)<1e-14,"profile native-grid scaling");
  require(unit.length_unit==3*std::numbers::pi/16,"observer length changed");
  const auto radius=s.R;s.projectHydrostatic();
  for(int f=0;f<NF;++f)require((s.R[f]-radius[f]).cwiseAbs().maxCoeff()<1e-8,"fitted profile not hydrostatic");
  // Verify the ACTUAL stock source, not the unavailable relative-U stencil.
  s.param.c1.fill(0);s.param.c2.fill(0);s.Deltat=.02;
  for(int f=0;f<NF;++f){s.Rho[f].setConstant(1+f);s.U[f].setConstant(.1+.03*f);}
  const auto old=s.U;s.solveConductionLAPACKE();
  for(int f=0;f<NF;++f)for(int i=1;i<199;++i) {
    double heating=0;
    for(int h=0;h<NF;++h) {
      if(f!=FB && h!=FB)continue;
      const int donor=f==FB?h:f;
      heating+=s.param.c4[f*NF+h]*s.Rho[h][i]*(1.5/std::sqrt(old[donor][i])-
        .5*s.U[donor][i]/std::pow(old[donor][i],1.5));
    }
    require(std::abs((s.U[f][i]-old[f][i])/s.Deltat-heating)<1e-12,"stock heating mismatch");
    require(s.U[f][199]==old[f][199],"stock outer U boundary changed");
  }
  bool threw=false;try{initializeHydroHeggie(s,p,{0,1});}catch(const std::invalid_argument&){threw=true;}
  require(threw,"invalid scale accepted");
  std::cout<<"stock Heggie initialization, native boundaries, heating and evolve checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
