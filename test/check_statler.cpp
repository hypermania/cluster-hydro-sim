#include "../src/reproduction.hpp"
#include "../src/statler_observer.hpp"
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool value,const char* text){if(!value)throw std::runtime_error(text);}
void near(double a,double b,double rel=2e-10){
  require(std::abs(a-b)<=rel*std::max({1e-15,std::abs(a),std::abs(b)}),"numeric regression");
}
struct Count {int calls=0;void operator()(const ThreeFluidSim&){++calls;}};

int main(){try {
  for(bool direct:{true,false}) {
    ThreeFluidSim s;auto p=statlerParameters(direct);p.maxSteps=100;
    initializeCaptureCluster(s,p);
    near(s.options.time_unit_over_trh,2.4284946,1e-7);
    StatlerObserver obs(s.parameters());s.evolve(obs);
    require(s.step==100 && obs.history.time_trh.size()==101,"accepted observer count");
    // Legacy 500-zone histories, recorded before moving the runner into evolve().
    near(obs.history.time_trh.back(),direct?0.15591076221655692:0.16209309361780602);
    near(obs.history.binary_number.back(),direct?3.807888406935049:3.9603263124390558);
    near(obs.history.central_density.back(),direct?0.25348735357006125:0.2540763025829066);
    near(obs.history.mass.back(),direct?0.3315486613243248:0.3315397526200759);
    near(obs.history.energy.back(),direct?-0.01617449405166721:-0.016173836939862837);
    const std::string dir=direct?"output/test_statler/direct/":"output/test_statler/control/";
    std::filesystem::create_directories(dir);s.saveParams(dir);obs.save(dir);
  }
  std::cout<<"legacy histories passed\n";
  ThreeFluidSim s;auto p=statlerParameters();p.N=150;p.maxSteps=1;
  initializeCaptureCluster(s,p);
  const auto oldrho=s.Rho,oldu=s.U;
  s.prepareFormationSource();const auto transfer=s.formationSource;
  s.applyBinaryFormation();
  for(int j=0;j<s.N;++j){
    near(s.Rho[FS][j]+s.Rho[FB][j],oldrho[FS][j]+oldrho[FB][j]);
    const double before=oldrho[FS][j]*oldu[FS][j]+oldrho[FB][j]*oldu[FB][j];
    const double after=s.Rho[FS][j]*s.U[FS][j]+s.Rho[FB][j]*s.U[FB][j];
    near(after,before-0.5*transfer[j]*oldu[FS][j]);
  }
  std::cout<<"formation invariants passed\n";
  initializeCaptureCluster(s,p);s.maxTime=0;Count zero;s.evolve(zero);
  require(s.step==0 && zero.calls==1,"equal endpoint must stop");
  initializeCaptureCluster(s,p);s.maxTime=1e-6;Count endpoint;s.evolve(endpoint);
  near(s.totalTime,1e-6);require(endpoint.calls==2,"endpoint observer count");
  std::cout<<"endpoint checks passed\n";
  initializeCaptureCluster(s,p);s.options.donor_fraction_limit=1e-9;
  Count retry;s.evolve(retry);
  require(s.rejected_steps>0 && retry.calls==2,"rejected trials leaked to observers");
  initializeCaptureCluster(s,p);s.options.donor_fraction_limit=1e-12;
  s.options.max_retries=0;
  const auto before_u=s.U,before_p=s.P;
  bool rejected=false;
  try{s.advanceAcceptedStep();}catch(const std::runtime_error&){rejected=true;}
  require(rejected && s.totalTime==0 && s.cumulative_formed_binaries==0,
          "failed trial changed accepted state");
  for(int f=0;f<NF;++f){
    require((s.U[f]-before_u[f]).norm()==0,"retry did not restore U");
    require((s.P[f]-before_p[f]).norm()==0,"retry did not restore P");
  }
  initializeCaptureCluster(s,p);s.binary_formation=0;s.options.timestep_controller=0;
  s.selectNextTimestep(0.001,0,s.centralDensity());
  require(std::isfinite(s.Deltat)&&s.Deltat>0,"zero error timestep");
  s.U[0][0]=std::numeric_limits<double>::quiet_NaN();
  bool failed=false;try{s.sanityCheck();}catch(const std::runtime_error&){failed=true;}
  require(failed,"NaN must throw");
  std::cout<<"statler_shared_evolution_passed=1\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
