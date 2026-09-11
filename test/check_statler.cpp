#include "../src/statler_reproduction.hpp"
#include "../src/observer.hpp"
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}
void near(double a,double b,double rel=2e-10){
  if(!std::isfinite(a)||!std::isfinite(b)||
     std::abs(a-b)>rel*std::max({1e-15,std::abs(a),std::abs(b)}))
    throw std::runtime_error("numeric regression: "+std::to_string(a)+" vs "+std::to_string(b));
}
struct Count {int calls=0;void operator()(const ThreeFluidSim&){++calls;}};

int main(){try {
  const StatlerInitParam initial;
  const auto observing=statlerObserverParameters(initial);
  near(observing.time_unit_over_trh,2.4284946,1e-7);
  for(bool direct:{true,false}) {
    ThreeFluidSim s;s.param=statlerParameters();s.param.maxSteps=100;
    initializeCaptureCluster(s,initial);
    near(s.param.capture_coefficient,191.69523955225444);
    if(!direct)s.param.c4.fill(0);
    const auto before=s.param;
    const std::string dir=direct?"output/test_statler/direct/":"output/test_statler/control/";
    std::filesystem::create_directories(dir+"initialization");
    std::filesystem::create_directories(dir+"observer");
    s.saveParams(dir);
    save_param_for_Mathematica(initial,dir+"initialization/");
    save_param_for_Mathematica(observing,dir+"observer/");
    StatlerObserver obs(observing);Count count;ObserverPack pack(obs,count);s.evolve(pack);
    require(s.step==100 && obs.history.time_trh.size()==101 && count.calls==101,
            "composed accepted observers");
    require(s.param.Deltat==before.Deltat && s.Deltat!=s.param.Deltat,
            "adaptive timestep mutated configuration");
    // First-step legacy values verify unchanged physical operators. Later
    // checkpoints use the requested U-only controller, not the removed clamp.
    near(obs.history.binary_number[1],direct?0.05921317884254938:0.05921317884250146);
    near(obs.history.energy[1],direct?-0.01633436448083242:-0.016334364480814084);
    near(obs.history.time_trh.back(),direct?0.15480564391023888:0.1614706543170104);
    near(obs.history.binary_number.back(),direct?3.7803196421511673:3.944598556969823);
    near(obs.history.central_density.back(),direct?0.25338476565548906:0.2540186390001849);
    near(obs.history.mass.back(),direct?0.3315133739660246:0.33150299314851595);
    near(obs.history.energy.back(),direct?-0.016171023627247684:-0.01617020194595756);
    obs.save(dir);
  }
  std::cout<<"legacy histories and independent parameters passed\n";
  ThreeFluidSim s;
  auto reset=[&]{s.param=statlerParameters();s.param.N=150;s.param.maxSteps=1;
                 initializeCaptureCluster(s,initial);};
  reset();
  const auto oldrho=s.Rho,oldu=s.U;
  s.applyBinaryFormation();
  for(int j=0;j<s.param.N;++j){
    const double transfer=s.param.mb*s.param.capture_coefficient*
      std::pow(oldrho[FS][j],2)/std::pow(oldu[FS][j],0.6)*s.Deltat;
    near(s.Rho[FS][j]+s.Rho[FB][j],oldrho[FS][j]+oldrho[FB][j]);
    const double before=oldrho[FS][j]*oldu[FS][j]+oldrho[FB][j]*oldu[FB][j];
    const double after=s.Rho[FS][j]*s.U[FS][j]+s.Rho[FB][j]*s.U[FB][j];
    near(after,before-0.5*transfer*oldu[FS][j]);
  }
  reset();s.param.maxTime=0;Count zero;s.evolve(zero);
  require(s.step==0 && zero.calls==1,"equal endpoint must stop");
  reset();s.param.maxTime=1e-6;Count endpoint;s.evolve(endpoint);
  near(s.totalTime,1e-6);require(endpoint.calls==2,"endpoint observer count");
  reset();s.param.capture_coefficient=1e30;
  Count failure;bool aborted=false;
  try{s.evolve(failure);}catch(const std::runtime_error&){aborted=true;}
  require(aborted && failure.calls==1 && s.step==0 && s.totalTime==0,
          "invalid formation did not terminate before accepting or observing the step");
  near(s.Deltat,s.param.Deltat); // No silent retry with a smaller step.
  reset();
  s.param.capture_coefficient*=0.01*s.Rho[FS][0]/(s.param.mb*s.captureNumberRate(0)*s.Deltat);
  const double donor=s.Rho[FS][0];s.applyBinaryFormation();
  near(s.Rho[FS][0],0.99*donor); // No residual 0.5% donor limit.
  reset();s.param.binary_formation=0;
  s.selectNextTimestep(0.001,0);near(s.Deltat,s.param.max_timestep);
  s.selectNextTimestep(0.001,1);near(s.Deltat,1e-6);
  s.Rho[0][0]*=100;
  s.selectNextTimestep(0.001,1);near(s.Deltat,1e-6); // U only.
  reset();s.previousU=s.U;s.U[0][0]=std::numeric_limits<double>::quiet_NaN();
  aborted=false;
  try{s.conductionChange();}catch(const std::runtime_error&){aborted=true;}
  require(aborted,"NaN energy change was hidden by the maximum reduction");
  reset();s.param.maxTime=0;s.U[0][0]=std::numeric_limits<double>::quiet_NaN();
  bool failed=false;Count validation;
  try{s.evolve(validation);}catch(const std::runtime_error&){failed=true;}
  require(failed && validation.calls==0,"enabled validation must reject NaN before observing");
  s.param.runtime_validation=0;s.evolve(validation);
  require(validation.calls==1,"disabled full-state validation still scanned the state");
  // Disabling scans does not disable structural controls or solver errors.
  s.param.max_timestep=0;failed=false;
  try{s.evolve(validation);}catch(const std::invalid_argument&){failed=true;}
  require(failed,"invalid controls accepted with scans disabled");
  reset();s.param.runtime_validation=0;++s.param.N;failed=false;
  try{s.evolve(validation);}catch(const std::invalid_argument&){failed=true;}
  require(failed,"grid configuration differs from allocated storage");
  reset();Count enabled;s.evolve(enabled);const auto checked=s.U;
  reset();s.param.runtime_validation=0;Count disabled;s.evolve(disabled);
  for(int f=0;f<NF;++f)require((checked[f]-s.U[f]).norm()==0,"validation changed numerics");
  std::cout<<"statler_shared_evolution_passed=1\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
