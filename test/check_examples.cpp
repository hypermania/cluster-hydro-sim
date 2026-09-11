// Bounded smoke checks for the eight existing main.cpp sample configurations.
#include "../src/observer.hpp"
#include <stdexcept>

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
    if(mode==6)p.c4.fill(0);
  }
  p.maxSteps=100;
}

int main(){try {
  for(int mode=0;mode<8;++mode) {
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
