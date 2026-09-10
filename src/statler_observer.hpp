#pragma once
#include "three_fluid.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// Reference-unit conversions and diagnostic approximations belong to observers,
// not the simulator's dimensionless internal state.
namespace statler_diagnostics {
struct Rates {
  double capture_number_per_hydro_time = 0.0;
  double three_body_number_per_hydro_time = 0.0;
  double capture_energy_per_hydro_time = 0.0;
  double direct_single_from_binary_single = 0.0;
  double direct_binary_from_binary_single = 0.0;
  double direct_binary_from_binary_binary = 0.0;
};

inline double shell_volume(const ThreeFluidSim &sim, const int fluid,
                    const int zone) {
  const double outer = sim.R[fluid][zone];
  const double inner = zone == 0 ? 0.0 : sim.R[fluid][zone - 1];
  return (std::pow(outer, 3) - std::pow(inner, 3)) / 3.0;
}

inline double total_mass(const ThreeFluidSim &sim) {
  double result = 0.0;
  for (int fluid = 0; fluid < NF; ++fluid) {
    result += sim.Menc[fluid][sim.N - 1];
  }
  return result;
}

inline double total_energy(const ThreeFluidSim &sim) {
  double random_energy = 0.0;
  double gravitational_energy = 0.0;
  double enclosed_mass = 0.0;
  for (int zone = 0; zone < sim.N; ++zone) {
    const double inner = zone == 0 ? 0.0 : sim.R[FS][zone - 1];
    const double outer = sim.R[FS][zone];
    double density = 0.0;
    for (int fluid = 0; fluid < NF; ++fluid) {
      const double cell_mass = sim.Rho[fluid][zone]
        * shell_volume(sim, fluid, zone);
      random_energy += cell_mass * sim.U[fluid][zone];
      density += sim.Rho[fluid][zone];
    }
    const double delta2 = outer * outer - inner * inner;
    const double delta3 = std::pow(outer, 3) - std::pow(inner, 3);
    const double delta5 = std::pow(outer, 5) - std::pow(inner, 5);
    gravitational_energy -= enclosed_mass * density * delta2 / 2.0;
    gravitational_energy -= density * density / 3.0
      * (delta5 / 5.0 - std::pow(inner, 3) * delta2 / 2.0);
    enclosed_mass += density * delta3 / 3.0;
  }
  return random_energy + gravitational_energy;
}

inline double half_mass_radius(const ThreeFluidSim &sim) {
  const double target = 0.5 * total_mass(sim);
  double inner_mass = 0.0;
  double inner_radius = 0.0;
  for (int zone = 0; zone < sim.N; ++zone) {
    double outer_mass = 0.0;
    for (int fluid = 0; fluid < NF; ++fluid) {
      outer_mass += sim.Menc[fluid][zone];
    }
    if (outer_mass >= target) {
      const double fraction = (target - inner_mass)
        / std::max(outer_mass - inner_mass,
                   std::numeric_limits<double>::min());
      const double outer_radius = sim.R[FS][zone];
      return std::cbrt(std::pow(inner_radius, 3)
                       + fraction * (std::pow(outer_radius, 3)
                                     - std::pow(inner_radius, 3)));
    }
    inner_mass = outer_mass;
    inner_radius = sim.R[FS][zone];
  }
  return sim.R[FS][sim.N - 1];
}

inline double core_radius(const ThreeFluidSim &sim) {
  double density = 0.0;
  double random_energy_density = 0.0;
  for (int fluid = 0; fluid < NF; ++fluid) {
    density += sim.Rho[fluid][0];
    random_energy_density += sim.Rho[fluid][0] * sim.U[fluid][0];
  }
  // Statler et al. equation (2.31): rc^2=3 v0^2/(4 pi G rho0).
  // In HydroSim units v0^2=2<U> and 4 pi G rho_unit r0^2=u_unit.
  return std::sqrt(6.0 * random_energy_density / (density * density));
}


inline double three_body_rate_density(const ThreeFluidSim &sim, const int zone) {
  const double reference_coulomb_log = sim.options.reference_coulomb_log;
  constexpr double code_mass_unit = 4.0 * std::numbers::pi;
  return 0.0009373511756007407 * (sim.ms * code_mass_unit)
    * std::pow(sim.Rho[FS][zone], 3)
    / (reference_coulomb_log * std::pow(sim.U[FS][zone], 4.5));
}

inline Rates evaluate_rates(const ThreeFluidSim &sim) {
  Rates rates;
  const double capture_coefficient = sim.options.capture_coefficient;
  for (int zone = 0; zone < sim.N; ++zone) {
    const double volume = shell_volume(sim, FS, zone);
    const double capture_number_density_rate = capture_coefficient
      * std::pow(sim.Rho[FS][zone], 2)
      / std::pow(sim.U[FS][zone], 0.6);
    rates.capture_number_per_hydro_time +=
      capture_number_density_rate * volume;
    rates.three_body_number_per_hydro_time +=
      three_body_rate_density(sim, zone) * volume;

    // A captured binary is inserted with U_b=U_s/2, so the resolved
    // translational random energy loses half of the transferred mass times U_s.
    rates.capture_energy_per_hydro_time -= 0.5 * sim.mb
      * capture_number_density_rate * sim.U[FS][zone] * volume;

    const double rho_s = sim.Rho[FS][zone];
    const double rho_b = sim.Rho[FB][zone];
    const double root_u_s = std::sqrt(sim.U[FS][zone]);
    const double root_u_b = std::sqrt(sim.U[FB][zone]);
    rates.direct_single_from_binary_single += rho_s * sim.c4[FS][FB]
      * rho_b / root_u_s * volume;
    rates.direct_binary_from_binary_single += rho_b * sim.c4[FB][FS]
      * rho_s / root_u_s * volume;
    rates.direct_binary_from_binary_binary += rho_b * sim.c4[FB][FB]
      * rho_b / root_u_b * volume;
  }
  return rates;
}


struct History {
  std::vector<double> time_trh;
  std::vector<double> central_density;
  std::vector<double> central_single_velocity_squared;
  std::vector<double> core_radius;
  std::vector<double> half_mass_radius;
  std::vector<double> binary_number;
  std::vector<double> central_density_ratio;
  std::vector<double> mean_binary_age_trh;
  std::vector<double> capture_number_rate_per_year;
  std::vector<double> three_body_number_rate_per_year;
  std::vector<double> capture_energy_rate;
  std::vector<double> direct_single_from_binary_single;
  std::vector<double> direct_binary_from_binary_single;
  std::vector<double> direct_binary_from_binary_binary;
  std::vector<double> mass;
  std::vector<double> energy;

  void record(const ThreeFluidSim &sim, const double hydro_time_per_trh,
              const double hydro_time_myr, const double binary_age_moment_trh) {
    const double current_time_trh = sim.totalTime * hydro_time_per_trh;
    const double rho0 = sim.Rho[FS][0] + sim.Rho[FB][0] + sim.Rho[FD][0];
    const double n_binary = sim.Menc[FB][sim.N - 1] / sim.mb;
    const Rates rates = evaluate_rates(sim);
    const double energy_rate_conversion = 9.0 / hydro_time_per_trh;

    time_trh.push_back(current_time_trh);
    central_density.push_back(3.0 * rho0 / (4.0 * std::numbers::pi));
    central_single_velocity_squared.push_back(6.0 * sim.U[FS][0]);
    core_radius.push_back(statler_diagnostics::core_radius(sim));
    half_mass_radius.push_back(statler_diagnostics::half_mass_radius(sim));
    binary_number.push_back(n_binary);
    central_density_ratio.push_back(sim.Rho[FB][0] / sim.Rho[FS][0]);
    mean_binary_age_trh.push_back(n_binary > 0.0
      ? binary_age_moment_trh / n_binary : 0.0);
    capture_number_rate_per_year.push_back(
      rates.capture_number_per_hydro_time / (hydro_time_myr * 1.0e6));
    three_body_number_rate_per_year.push_back(
      rates.three_body_number_per_hydro_time / (hydro_time_myr * 1.0e6));
    capture_energy_rate.push_back(
      rates.capture_energy_per_hydro_time * energy_rate_conversion);
    direct_single_from_binary_single.push_back(
      rates.direct_single_from_binary_single * energy_rate_conversion);
    direct_binary_from_binary_single.push_back(
      rates.direct_binary_from_binary_single * energy_rate_conversion);
    direct_binary_from_binary_binary.push_back(
      rates.direct_binary_from_binary_binary * energy_rate_conversion);
    mass.push_back(total_mass(sim));
    energy.push_back(total_energy(sim));
  }

  void save(const std::string &directory) const {
    write_to_file(time_trh, directory + "history_time_trh.dat");
    write_to_file(central_density, directory + "history_rho0.dat");
    write_to_file(central_single_velocity_squared,
                  directory + "history_vms2_0.dat");
    write_to_file(core_radius, directory + "history_rc.dat");
    write_to_file(half_mass_radius, directory + "history_rh.dat");
    write_to_file(binary_number, directory + "history_nb.dat");
    write_to_file(central_density_ratio, directory + "history_ratio.dat");
    write_to_file(mean_binary_age_trh,
                  directory + "history_mean_binary_age_trh.dat");
    write_to_file(capture_number_rate_per_year,
                  directory + "history_capture_number_rate_per_year.dat");
    write_to_file(three_body_number_rate_per_year,
                  directory + "history_threebody_number_rate_per_year.dat");
    write_to_file(capture_energy_rate,
                  directory + "history_capture_energy_rate.dat");
    write_to_file(direct_single_from_binary_single,
                  directory + "history_direct_s_bs.dat");
    write_to_file(direct_binary_from_binary_single,
                  directory + "history_direct_b_bs.dat");
    write_to_file(direct_binary_from_binary_binary,
                  directory + "history_direct_b_bb.dat");
    write_to_file(mass, directory + "history_mass.dat");
    write_to_file(energy, directory + "history_energy.dat");
  }
};

struct Snapshot {
  double time_trh = 0.0;
  std::vector<double> radius;
  std::vector<double> rho_single;
  std::vector<double> rho_binary;
  std::vector<double> u_single;
  std::vector<double> u_binary;
  std::vector<double> luminosity;
};

inline Snapshot make_snapshot(const ThreeFluidSim &sim, const double time_trh,
                       const double hydro_time_per_trh) {
  Snapshot snapshot;
  snapshot.time_trh = time_trh;
  snapshot.radius.resize(sim.N);
  snapshot.rho_single.resize(sim.N);
  snapshot.rho_binary.resize(sim.N);
  snapshot.u_single.resize(sim.N);
  snapshot.u_binary.resize(sim.N);
  snapshot.luminosity.assign(sim.N, 0.0);
  const double luminosity_conversion = 9.0 / hydro_time_per_trh;
  for (int zone = 0; zone < sim.N; ++zone) {
    snapshot.radius[zone] = sim.R[FS][zone];
    snapshot.rho_single[zone] = sim.Rho[FS][zone];
    snapshot.rho_binary[zone] = sim.Rho[FB][zone];
    snapshot.u_single[zone] = sim.U[FS][zone];
    snapshot.u_binary[zone] = sim.U[FB][zone];
    if (zone == 0) {
      continue;
    }
    for (int fluid = 0; fluid < NF; ++fluid) {
      snapshot.luminosity[zone] -= sim.c2[fluid] * sim.Rho[fluid][zone]
        * std::pow(sim.R[fluid][zone], 2)
        * (std::sqrt(sim.U[fluid][zone])
           - std::sqrt(sim.U[fluid][zone - 1]))
        / (sim.R[fluid][zone] - sim.R[fluid][zone - 1]);
    }
    snapshot.luminosity[zone] *= luminosity_conversion;
  }
  return snapshot;
}

struct SnapshotStore {
  std::vector<double> targets;
  std::vector<Snapshot> snapshots;
  std::size_t next_target = 0;
  Snapshot peak;
  double peak_density = -1.0;
  bool save_peak = true;

  void observe(const ThreeFluidSim &sim, const double hydro_time_per_trh) {
    const double time_trh = sim.totalTime * hydro_time_per_trh;
    while (next_target < targets.size() && time_trh >= targets[next_target]) {
      snapshots.push_back(make_snapshot(sim, time_trh, hydro_time_per_trh));
      ++next_target;
    }
    const double central_density = sim.Rho[FS][0] + sim.Rho[FB][0];
    if (central_density > peak_density) {
      peak_density = central_density;
      peak = make_snapshot(sim, time_trh, hydro_time_per_trh);
    }
  }

  void save(const std::string &directory, const int zones) const {
    std::vector<Snapshot> output = snapshots;
    if (save_peak && peak_density >= 0) output.push_back(peak);
    std::sort(output.begin(), output.end(), [](const Snapshot &left,
                                                const Snapshot &right) {
      return left.time_trh < right.time_trh;
    });
    std::vector<double> times;
    std::vector<double> radius;
    std::vector<double> rho_single;
    std::vector<double> rho_binary;
    std::vector<double> u_single;
    std::vector<double> u_binary;
    std::vector<double> luminosity;
    for (const Snapshot &snapshot : output) {
      times.push_back(snapshot.time_trh);
      radius.insert(radius.end(), snapshot.radius.begin(), snapshot.radius.end());
      rho_single.insert(rho_single.end(), snapshot.rho_single.begin(),
                        snapshot.rho_single.end());
      rho_binary.insert(rho_binary.end(), snapshot.rho_binary.begin(),
                        snapshot.rho_binary.end());
      u_single.insert(u_single.end(), snapshot.u_single.begin(),
                      snapshot.u_single.end());
      u_binary.insert(u_binary.end(), snapshot.u_binary.begin(),
                      snapshot.u_binary.end());
      luminosity.insert(luminosity.end(), snapshot.luminosity.begin(),
                        snapshot.luminosity.end());
    }
    write_to_file(times, directory + "snapshot_time_trh.dat");
    write_to_file(radius, directory + "snapshot_radius.dat");
    write_to_file(rho_single, directory + "snapshot_rho_s.dat");
    write_to_file(rho_binary, directory + "snapshot_rho_b.dat");
    write_to_file(u_single, directory + "snapshot_u_s.dat");
    write_to_file(u_binary, directory + "snapshot_u_b.dat");
    write_to_file(luminosity, directory + "snapshot_luminosity.dat");
    write_to_file(std::vector<long long>{zones},
                  directory + "snapshot_zone_count.dat");
  }
};


} // namespace statler_diagnostics

struct StatlerObserver {
  statler_diagnostics::History history;
  statler_diagnostics::SnapshotStore snapshots;
  ThreeFluidParam config;
  double age_moment_hydro = 0;
  double previous_time = 0, previous_number = 0, previous_formed = 0;
  bool initialized = false;
  explicit StatlerObserver(const ThreeFluidParam& p):config(p) {
    if(!p.statler_observer || p.history_stride<1 || p.output_format!=1 ||
       p.snapshot_count<0 || p.snapshot_count>16 ||
       !std::isfinite(p.time_unit_myr) || p.time_unit_myr<=0 ||
       !std::isfinite(p.time_unit_over_trh) || p.time_unit_over_trh<=0)
      throw std::invalid_argument("invalid Statler observer configuration");
    snapshots.targets.assign(p.snapshot_times_trh.begin(),
                              p.snapshot_times_trh.begin()+p.snapshot_count);
    for(double t:snapshots.targets)
      if(!std::isfinite(t)||t<0) throw std::invalid_argument("invalid snapshot time");
    if(!std::is_sorted(snapshots.targets.begin(),snapshots.targets.end()))
      throw std::invalid_argument("unordered snapshot targets");
    snapshots.save_peak=p.save_peak_snapshot!=0;
  }
  void operator()(const ThreeFluidSim& sim) {
    const double number=sim.Menc[FB][sim.N-1]/sim.mb;
    if(initialized) {
      const double dt=sim.totalTime-previous_time;
      const double formed=sim.cumulative_formed_binaries-previous_formed;
      // Aggregate age estimate; numerical remapping loss is assumed age-neutral.
      age_moment_hydro+=previous_number*dt+0.5*formed*dt;
      if(previous_number+formed>0) age_moment_hydro*=number/(previous_number+formed);
    }
    previous_time=sim.totalTime;previous_number=number;
    previous_formed=sim.cumulative_formed_binaries;initialized=true;
    if(sim.step%config.history_stride==0 || sim.stopCondition())
      history.record(sim,config.time_unit_over_trh,config.time_unit_myr,
                     age_moment_hydro*config.time_unit_over_trh);
    snapshots.observe(sim,config.time_unit_over_trh);
  }
  void save(const std::string& directory) const {
    history.save(directory);
    snapshots.save(directory,static_cast<int>(config.N));
    write_to_file(std::vector<double>{config.time_unit_myr},directory+"hydro_time_myr.dat");
    write_to_file(std::vector<double>{config.time_unit_over_trh},directory+"hydro_time_per_trh.dat");
  }
};
