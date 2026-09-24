#pragma once
#include "moving_three_fluid.hpp"
#include "statler_observer.hpp"

// Read-only diagnostic adapter. The scratch ThreeFluidSim is never evolved:
// it supplies the established Statler units/schema without a second integrator.
struct MovingStatlerObserver {
  StatlerObserver diagnostic;
  ThreeFluidSim view;
  std::ofstream budget;
  explicit MovingStatlerObserver(const StatlerObserverParam& p,const std::string& dir)
      :diagnostic(p),budget(dir+"budget.csv") {
    budget.exceptions(std::ios::badbit|std::ios::failbit);
    budget<<std::setprecision(17)<<"step,time,dt,formed,capture_energy,mass_plus_outflow,energy_ledger_error,last_change\n";
  }
  void profile(const MovingThreeFluidSim& s,statler_diagnostics::Snapshot& p) {
    auto luminosity=[&](int face) {
      if(face==0)return 0.;
      double result=0;
      for(int f=0;f<NF;++f) {
        if(face==s.zones()) {
          const int i=face-1;
          result+=s.param.c2[f]*std::pow(s.faces()[face],2)*s.value(i,f,s.RHO)*
            std::sqrt(s.value(i,f,s.U))/(s.param.thermal_length_over_radius*s.faces()[face]+
                                       s.faces()[face]-s.radius(i));
        } else {
          const double l=s.value(face-1,f,s.RHO),r=s.value(face,f,s.RHO);
          result-=s.param.c2[f]*std::pow(s.faces()[face],2)*2*l*r/(l+r)*
            (std::sqrt(s.value(face,f,s.U))-std::sqrt(s.value(face-1,f,s.U)))/
            (s.radius(face)-s.radius(face-1));
        }
      }
      return result*9/diagnostic.param.time_unit_over_trh;
    };
    for(int i=0;i<s.zones();++i) {
      p.radius[i]=s.radius(i);
      // Snapshot estimate interpolated to the density/temperature centre.
      // This is not the accepted affine face flux used in the energy ledger.
      p.luminosity[i]=(luminosity(i)+luminosity(i+1))/2;
    }
  }
  void operator()(const MovingThreeFluidSim& s) {
    if(view.R[FS].size()!=s.zones())view.initSolver(s.zones());
    view.param.ms=s.param.mass[FS];view.param.mb=s.param.mass[FB];view.param.md=s.param.mass[FD];
    view.param.c1=s.param.c1;view.param.c2=s.param.c2;view.param.c4=s.param.c4;
    view.param.capture_coefficient=s.param.capture_coefficient;
    view.totalTime=s.totalTime;view.step=s.step;
    view.cumulative_formed_binaries=s.cumulative_formed_binaries;
    double bulk_energy=0,mass=0;
    for(int f=0;f<NF;++f) {
      mass+=s.value(s.zones()-1,f,s.MASS)+s.escaped_mass[f];
      for(int i=0;i<s.zones();++i) {
        view.R[f][i]=s.faces()[i+1];view.Rho[f][i]=s.value(i,f,s.RHO);
        view.U[f][i]=s.value(i,f,s.U);view.Menc[f][i]=s.value(i,f,s.MASS);
        bulk_energy+=0.5*s.param.epsilon*s.volume(i)*s.value(i,f,s.RHO)*std::pow(s.value(i,f,s.VEL),2);
      }
    }
    const auto count=diagnostic.snapshots.snapshots.size();
    diagnostic(view);
    diagnostic.history.energy.back()+=bulk_energy;
    for(auto j=count;j<diagnostic.snapshots.snapshots.size();++j)
      profile(s,diagnostic.snapshots.snapshots[j]);
    if(diagnostic.snapshots.peak.time_trh==s.totalTime*diagnostic.param.time_unit_over_trh)
      profile(s,diagnostic.snapshots.peak);
    budget<<s.step<<','<<s.totalTime<<','<<s.Deltat<<','<<s.cumulative_formed_binaries
      <<','<<s.cumulative_capture_energy<<','<<mass<<','<<s.energy_ledger_error<<','<<s.last_change<<'\n';
    if(s.step%10000==0) {
      budget.flush();
      std::cout<<"moving Statler step="<<s.step<<" time_trh="
        <<s.totalTime*diagnostic.param.time_unit_over_trh<<" dt="<<s.Deltat
        <<" rho_s="<<s.value(0,FS,s.RHO)<<std::endl;
    }
  }
  void save(const std::string& dir) {diagnostic.save(dir);budget.flush();}
};
