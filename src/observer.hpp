#pragma once

#include "three_fluid.hpp"

struct KeyValueObserver {
  std::vector<double> t_list;
  std::array<std::vector<double>, NF> central_Rho;
  std::array<std::vector<double>, NF> central_U;
  std::array<std::vector<double>, NF> total_M;

  KeyValueObserver() {}
  
  void operator()(const ThreeFluidSim &sim) {
    const auto &Rho = sim.Rho;
    const auto &U = sim.U;
    const auto &Menc = sim.Menc;
    const double t = sim.totalTime;
    t_list.push_back(t);
    for(int f = 0; f < NF; ++f){
      central_Rho[f].push_back(Rho[f][0]);
      central_U[f].push_back(U[f][0]);
      total_M[f].push_back(Menc[f][sim.param.N-1]);
    }
  }

  void save(const std::string &dir) const {
    write_to_file(t_list, dir + "central_t.dat");
    write_to_file(central_Rho[FS], dir + "central_Rho_s.dat");
    write_to_file(central_Rho[FB], dir + "central_Rho_b.dat");
    write_to_file(central_Rho[FD], dir + "central_Rho_d.dat");
    
    write_to_file(central_U[FS], dir + "central_U_s.dat");
    write_to_file(central_U[FB], dir + "central_U_b.dat");
    write_to_file(central_U[FD], dir + "central_U_d.dat");
    
    write_to_file(total_M[FS], dir + "total_M_s.dat");
    write_to_file(total_M[FB], dir + "total_M_b.dat");
    write_to_file(total_M[FD], dir + "total_M_d.dat");
  }
};



struct ApproximateTimeObserver {
  std::vector<double> times;
  
  std::vector<double> t_list;
  std::vector<double> R_packed;
  std::vector<double> Rho_s_packed;
  std::vector<double> Rho_b_packed;
  std::vector<double> Rho_d_packed;
  std::vector<double> U_s_packed;
  std::vector<double> U_b_packed;
  std::vector<double> U_d_packed;

  ApproximateTimeObserver(const std::vector<double> &times_) : times(times_) {}
  
  void operator()(const ThreeFluidSim &sim) {
    const Eigen::VectorXd &R = sim.R[FS];
    const auto &Rho = sim.Rho;
    const auto &U = sim.U;
    const double t = sim.totalTime;
    const size_t current_idx = t_list.size();
    if(current_idx < times.size() && t >= times[current_idx]) {
      t_list.push_back(t);
      R_packed.insert(R_packed.end(), R.begin(), R.end());
      Rho_s_packed.insert(Rho_s_packed.end(), Rho[FS].begin(), Rho[FS].end());
      Rho_b_packed.insert(Rho_b_packed.end(), Rho[FB].begin(), Rho[FB].end());
      Rho_d_packed.insert(Rho_d_packed.end(), Rho[FD].begin(), Rho[FD].end());
      U_s_packed.insert(U_s_packed.end(), U[FS].begin(), U[FS].end());
      U_b_packed.insert(U_b_packed.end(), U[FB].begin(), U[FB].end());
      U_d_packed.insert(U_d_packed.end(), U[FD].begin(), U[FD].end());
    }
  }
  
  void save(const std::string &dir) const {
    write_to_file(t_list, dir + "t.dat");
    write_to_file(R_packed, dir + "R.dat");
    write_to_file(Rho_s_packed, dir + "Rho_s.dat");
    write_to_file(Rho_b_packed, dir + "Rho_b.dat");
    write_to_file(Rho_d_packed, dir + "Rho_d.dat");
    write_to_file(U_s_packed, dir + "U_s.dat");
    write_to_file(U_b_packed, dir + "U_b.dat");
    write_to_file(U_d_packed, dir + "U_d.dat");
  }
};


struct LagrangianRadiiObserver {
  std::vector<double> fraction;
  
  std::vector<double> t_list;
  std::array<std::vector<double>, NF> radii_packed; // fraction.size * t_list.size radii for each elem

  LagrangianRadiiObserver(const std::vector<double> &fraction_) : fraction(fraction_) {}
  
  void operator()(const ThreeFluidSim &sim) {
    const Eigen::VectorXd &R = sim.R[FS];
    const auto &Rho = sim.Rho;
    const auto &Menc = sim.Menc;
    const double t = sim.totalTime;
    const long long int N = sim.param.N;
    
    //const size_t current_idx = t_list.size();
    t_list.push_back(t);
    for(int f = 0; f < NF; ++f){
      for(double frac : fraction){
	double threshold = Menc[f][N-1] * frac;
	int idx = std::lower_bound(Menc[f].begin(), Menc[f].end(), threshold) - Menc[f].begin();
	double r = 0;
	if(idx == 0){
	  // pow(R[f][0], 3) / 3.0 * Rho[f][0] == threshold
	  r = pow(3.0 * threshold /  Rho[f][0], 1.0 / 3.0);
	} else if(idx < N) {
	  // Assume constant density in Lagrangian zone
	  // threshold - Menc[idx-1] == (pow(r, 3) - pow(R[f][idx-1], 3)) / 3.0 * Rho[f][idx]
	  r = pow(3.0 / Rho[f][idx] * (threshold - Menc[f][idx-1]) + pow(R[idx-1], 3), 1.0 / 3.0);
	  
	  // double t = (log(R[f][idx]) - log(R[f][idx-1])) / (log(R[f][idx]) - log(R[f][idx-1]));
	  // r = exp(log(Rho[f][idx-1]) * (1-t) + log(Rho[f][idx]) * t);
	} else {
	  r = R[N-1];
	}
	radii_packed[f].push_back(r);
      }
    }
  }
  
  void save(const std::string &dir) const {
    write_to_file(fraction, dir + "radii_fraction.dat");
    write_to_file(t_list, dir + "radii_t.dat");
    write_to_file(radii_packed[FS], dir + "radii_s.dat");
    write_to_file(radii_packed[FB], dir + "radii_b.dat");
    write_to_file(radii_packed[FD], dir + "radii_d.dat");
  }
};


template<typename... Observers>
struct ObserverPack {
  std::tuple<Observers & ...> observers;
  
  ObserverPack(Observers & ... observers_) : observers(observers_...) {}

  void operator()(const ThreeFluidSim &sim) {
    std::apply([&](auto &&... args) { ((args(sim)), ...); }, observers);
  }

  void save(const std::string &dir) const {
    // std::apply([&](auto &&... args) { ((args.save(dir)), ...); }, observers);
    std::apply([&](auto &&... args) {
      ([&](auto &&arg){
	if constexpr (requires { arg.save(dir); }) {
	  arg.save(dir);
	}
      }(args), ...);
    }, observers);
  }
};
