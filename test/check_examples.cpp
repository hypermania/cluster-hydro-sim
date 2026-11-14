// Bounded sample checks, including the three-body-only formation mode.
#include "../src/observer.hpp"
#include <numbers>
#include <stdexcept>

void checkFormationModes() {
  ThreeFluidSim initial;initial.initSolver(8);
  initial.param.ms=.1;initial.param.mb=.2;initial.param.md=.01;
  initial.Deltat=.01;
  for(int f=0;f<NF;++f) {
    initial.R[f].setLinSpaced(8,.1,2.);
    initial.Rho[f].setConstant(f==FS?2.:.3);
    initial.U[f].setConstant(f==FS?.5:.2);
    initial.P[f]=(2./3.)*initial.Rho[f].array()*initial.U[f].array();
  }
  initial.updateEnclosedMass();
  const auto& p=initial.param;
  double mass=0;for(int f=0;f<NF;++f)mass+=initial.Menc[f][p.N-1];
  const double log_lambda=std::log(GAMMA_LAMBDA*2*mass/(p.ms+p.md));
  const double rate3=0.0009373511756007407*(p.ms*4*std::numbers::pi)*8/
                    (log_lambda*std::pow(.5,4.5));
  const double capture=3.52e-6*std::pow(p.ms*4*std::numbers::pi,-2)*4/
                       (log_lambda*std::sqrt(.5));
  auto close=[](double actual,double expected) {
    if(!std::isfinite(actual) || std::abs(actual-expected)>2e-15+3e-12*std::abs(expected))
      throw std::runtime_error("formation rate or conserved state mismatch");
  };
  for(auto mode:{BINARY_FORMATION_MODE_2,BINARY_FORMATION_RODRIGUEZ}) {
    auto s=initial;s.param.binary_formation=mode;s.applyBinaryFormation();
    const double added=p.mb*initial.Deltat*(rate3+(mode==BINARY_FORMATION_MODE_2?capture:0));
    for(int j=0;j<p.N;++j) {
      close(s.Rho[FB][j]-.3,added);close(2-s.Rho[FS][j],added);
      close(s.U[FB][j],(.3*.2+added*.5)/(.3+added));
      close(s.Rho[FS][j]*s.U[FS][j]+s.Rho[FB][j]*s.U[FB][j],1.06);
      close(s.U[FS][j],.5);close(s.Rho[FD][j],.3);close(s.U[FD][j],.2);
      close(s.Menc[FS][j]+s.Menc[FB][j],initial.Menc[FS][j]+initial.Menc[FB][j]);
      for(int f=0;f<NF;++f) {
        close(s.P[f][j],(2./3.)*s.Rho[f][j]*s.U[f][j]);
        close(s.R[f][j],initial.R[f][j]);
      }
    }
  }
  std::cout<<"MODE_2 / RODRIGUEZ source and local conservation checks passed\n";
}

void setup(ThreeFluidSim& s,int mode) {
  auto& p=s.param;
  s.initSolver(mode<5?150:500);
  if(mode<5) {
    s.initCoeffsYiming();
    p.md=mode<3?p.ms:0.1*p.ms;
    if(mode<3)p.c2[FD]=p.c2[FS];
    s.initPlummerYiming(mode==0||mode==2?1.0:0.5,1e-10,
                        mode==0||mode==2?1e-10:1.0,1.0,1.0);
    if(mode==2||mode==4) {
      p.tidal_cutoff=1;p.tidal_radius=mode==2?10:2;
      p.tidal_cutoff_factor=mode==2?50:10;
    }
  } else {
    s.initPlummer(1,1e-10,1e-10,1,1);
    double mass=0;for(int f=0;f<NF;++f)mass+=s.Menc[f][p.N-1];
    p.ms=mass/1e6;p.mb=2*p.ms;p.md=1e-10*p.ms;
    s.initCoeffs(1e6,2,1e-10);
    if(mode>5)p.binary_formation=BINARY_FORMATION_MODE_2;
    if(mode==8)p.binary_formation=BINARY_FORMATION_RODRIGUEZ;
    if(mode==6)p.c4.fill(0);
  }
  p.maxSteps=100;
}

int main(){try {
  checkFormationModes();
  for(int mode=0;mode<9;++mode) {
    ThreeFluidSim s;setup(s,mode);
    KeyValueObserver observer;s.evolve(observer);
    if(s.step!=100 || observer.t_list.size()!=101)
      throw std::runtime_error("unexpected accepted history length");
    s.sanityCheck();
    for(int f=0;f<NF;++f)
      if(observer.total_M[f].back()!=s.Menc[f][s.param.N-1])
        throw std::runtime_error("observer mass mismatch");
    std::cout<<"example="<<mode<<" steps="<<s.step<<" time="<<s.totalTime<<" passed=1\n";
  }
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
