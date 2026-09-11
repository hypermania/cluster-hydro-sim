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
    // Legacy 500-zone histories; observers retain the reference conventions.
    near(obs.history.time_trh.back(),direct?0.15591076221655692:0.16209309361780602);
    near(obs.history.binary_number.back(),direct?3.807888406935049:3.9603263124390558);
    near(obs.history.central_density.back(),direct?0.25348735357006125:0.2540763025829066);
    near(obs.history.mass.back(),direct?0.3315486613243248:0.3315397526200759);
    near(obs.history.energy.back(),direct?-0.01617449405166721:-0.016173836939862837);
    obs.save(dir);
  }
  std::cout<<"legacy histories and independent parameters passed\n";
  ThreeFluidSim s;
  auto reset=[&]{s.param=statlerParameters();s.param.N=150;s.param.maxSteps=1;
                 initializeCaptureCluster(s,initial);};
  reset();
  const auto oldrho=s.Rho,oldu=s.U;
  s.prepareFormationSource();const auto transfer=s.formationSource;
  s.applyBinaryFormation();
  for(int j=0;j<s.param.N;++j){
    near(s.Rho[FS][j]+s.Rho[FB][j],oldrho[FS][j]+oldrho[FB][j]);
    const double before=oldrho[FS][j]*oldu[FS][j]+oldrho[FB][j]*oldu[FB][j];
    const double after=s.Rho[FS][j]*s.U[FS][j]+s.Rho[FB][j]*s.U[FB][j];
    near(after,before-0.5*transfer[j]*oldu[FS][j]);
  }
  reset();s.param.maxTime=0;Count zero;s.evolve(zero);
  require(s.step==0 && zero.calls==1,"equal endpoint must stop");
  reset();s.param.maxTime=1e-6;Count endpoint;s.evolve(endpoint);
  near(s.totalTime,1e-6);require(endpoint.calls==2,"endpoint observer count");
  reset();s.param.donor_fraction_limit=1e-9;
  Count retry;s.evolve(retry);
  require(s.rejected_steps>0 && retry.calls==2,"rejected trials leaked to observers");
  const double accepted_dt=s.totalTime;
  const auto retry_u=s.U,retry_p=s.P,retry_rho=s.Rho;
  reset();s.param.Deltat=accepted_dt;Count fresh;s.evolve(fresh);
  for(int f=0;f<NF;++f){
    require((s.U[f]-retry_u[f]).norm()==0,"retry U differs from fresh reduced step");
    require((s.P[f]-retry_p[f]).norm()==0,"retry P differs from fresh reduced step");
    require((s.Rho[f]-retry_rho[f]).norm()==0,"retry density differs from fresh reduced step");
  }
  reset();s.param.binary_formation=0;
  s.selectNextTimestep(0.001,0,s.centralDensity());near(s.Deltat,0.0015);
  s.selectNextTimestep(0.001,1,s.centralDensity());near(s.Deltat,0.0005);
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
