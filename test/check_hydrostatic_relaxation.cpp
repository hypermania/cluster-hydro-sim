#include "../src/three_fluid.hpp"

// Run independently from the repository root with:
//   make check_hydrostatic_relaxation
//   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./check_hydrostatic_relaxation

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

constexpr double residual_tolerance = 1e-10;
constexpr double modification_tolerance = 1e-12;
constexpr double invariant_tolerance = 1e-12;

bool finite_bits(const double value) {
  constexpr std::uint64_t exponent = UINT64_C(0x7ff0000000000000);
  return (std::bit_cast<std::uint64_t>(value) & exponent) != exponent;
}

double total_mass_at(const ThreeFluidSim &sim, const int i) {
  return sim.Menc[FS][i] + sim.Menc[FB][i] + sim.Menc[FD][i];
}

double shell_volume(const Eigen::VectorXd &radius, const int i) {
  const double inner_radius = (i == 0) ? 0.0 : radius[i-1];
  return (std::pow(radius[i], 3) - std::pow(inner_radius, 3)) / 3.0;
}

double relative_difference(const double value, const double reference) {
  return std::abs(value - reference)
    / std::max(std::abs(reference), std::numeric_limits<double>::min());
}

double hydrostatic_residual(const ThreeFluidSim &sim, const int fluid) {
  double maximum = 0.0;
  for (int i = 0; i < sim.param.N - 1; ++i) {
    const double span = (i == 0)
      ? sim.R[fluid][1]
      : sim.R[fluid][i+1] - sim.R[fluid][i-1];
    const double pressure_term = 4.0 * std::pow(sim.R[fluid][i], 2)
      * (sim.P[fluid][i+1] - sim.P[fluid][i]);
    const double gravity_term = total_mass_at(sim, i)
      * (sim.Rho[fluid][i] + sim.Rho[fluid][i+1]) * span;
    const double scale = std::abs(pressure_term) + std::abs(gravity_term);
    maximum = std::max(maximum,
      (scale > 0.0)
        ? std::abs(pressure_term + gravity_term) / scale
        : std::abs(pressure_term + gravity_term));
  }

  const int last = sim.param.N - 1;
  const double width = sim.R[fluid][last] - sim.R[fluid][last-1];
  const double pressure_term = -sim.P[fluid][last] / width;
  const double gravity_term = total_mass_at(sim, last)
    * sim.Rho[fluid][last] / std::pow(sim.R[fluid][last], 2);
  const double scale = std::abs(pressure_term) + std::abs(gravity_term);
  return std::max(maximum,
    (scale > 0.0)
      ? std::abs(pressure_term + gravity_term) / scale
      : std::abs(pressure_term + gravity_term));
}

void initialize_discrete_equilibrium(ThreeFluidSim &sim, const int zones) {
  sim.initSolver(zones);
  constexpr std::array<double, NF> central_density{0.5, 0.2, 0.3};
  constexpr std::array<double, NF> scale_radius{1.0, 0.7, 1.4};

  for (int f = 0; f < NF; ++f) {
    sim.R[f].array() = Eigen::pow(
      10.0, Eigen::VectorXd::LinSpaced(zones, -2.0, 2.0).array());
    for (int i = 0; i < zones; ++i) {
      const double shell_center = (i == 0)
        ? sim.R[f][0] / 2.0
        : (sim.R[f][i-1] + sim.R[f][i]) / 2.0;
      sim.Rho[f][i] = central_density[f]
        / std::pow(1.0 + std::pow(shell_center / scale_radius[f], 2), 2.5);
    }
  }
  sim.updateEnclosedMass();

  const int last = zones - 1;
  for (int f = 0; f < NF; ++f) {
    const double width = sim.R[f][last] - sim.R[f][last-1];
    sim.P[f][last] = total_mass_at(sim, last) * sim.Rho[f][last]
      * width / std::pow(sim.R[f][last], 2);
  }
  for (int i = zones - 2; i >= 0; --i) {
    for (int f = 0; f < NF; ++f) {
      const double span = (i == 0)
        ? sim.R[f][1]
        : sim.R[f][i+1] - sim.R[f][i-1];
      sim.P[f][i] = sim.P[f][i+1]
        + total_mass_at(sim, i)
          * (sim.Rho[f][i] + sim.Rho[f][i+1]) * span
          / (4.0 * std::pow(sim.R[f][i], 2));
    }
  }
  for (int f = 0; f < NF; ++f) {
    sim.U[f] = 1.5 * (sim.P[f].array() / sim.Rho[f].array()).matrix();
  }
}

bool check_grid(const int zones) {
  ThreeFluidSim sim;
  initialize_discrete_equilibrium(sim, zones);

  const auto initial_radius = sim.R;
  const auto initial_density = sim.Rho;
  const auto initial_pressure = sim.P;
  const auto initial_energy = sim.U;
  double initial_residual = 0.0;
  for (int f = 0; f < NF; ++f) {
    initial_residual = std::max(initial_residual,
                                hydrostatic_residual(sim, f));
  }

  for (int f = 0; f < NF; ++f) {
    sim.solveRelaxationLAPACKE(f);
    sim.solveRelaxationLAPACKE(f);
  }
  sim.updateEnclosedMass();

  double final_residual = 0.0;
  double maximum_radius_correction = 0.0;
  double maximum_state_change = 0.0;
  double maximum_mass_error = 0.0;
  double maximum_entropy_error = 0.0;
  bool finite_and_ordered = true;

  for (int f = 0; f < NF; ++f) {
    final_residual = std::max(final_residual,
                              hydrostatic_residual(sim, f));
    for (int i = 0; i < zones; ++i) {
      const double inner_width = (i == 0)
        ? initial_radius[f][0]
        : initial_radius[f][i] - initial_radius[f][i-1];
      const double outer_width = (i + 1 == zones)
        ? inner_width
        : initial_radius[f][i+1] - initial_radius[f][i];
      maximum_radius_correction = std::max(
        maximum_radius_correction,
        std::abs(sim.R[f][i] - initial_radius[f][i])
          / std::min(inner_width, outer_width));
      maximum_state_change = std::max({
        maximum_state_change,
        relative_difference(sim.Rho[f][i], initial_density[f][i]),
        relative_difference(sim.P[f][i], initial_pressure[f][i]),
        relative_difference(sim.U[f][i], initial_energy[f][i])});

      const double initial_mass = initial_density[f][i]
        * shell_volume(initial_radius[f], i);
      const double final_mass = sim.Rho[f][i] * shell_volume(sim.R[f], i);
      maximum_mass_error = std::max(
        maximum_mass_error, relative_difference(final_mass, initial_mass));
      const double initial_entropy = initial_pressure[f][i]
        / std::pow(initial_density[f][i], 5.0 / 3.0);
      const double final_entropy = sim.P[f][i]
        / std::pow(sim.Rho[f][i], 5.0 / 3.0);
      maximum_entropy_error = std::max(
        maximum_entropy_error,
        relative_difference(final_entropy, initial_entropy));

      finite_and_ordered = finite_and_ordered
        && finite_bits(sim.R[f][i]) && finite_bits(sim.Rho[f][i])
        && finite_bits(sim.P[f][i]) && finite_bits(sim.U[f][i])
        && sim.R[f][i] > 0.0 && sim.Rho[f][i] > 0.0
        && sim.P[f][i] > 0.0
        && (i == 0 || sim.R[f][i] > sim.R[f][i-1]);
    }
  }

  const bool passed = finite_and_ordered
    && initial_residual < residual_tolerance
    && final_residual < residual_tolerance
    && maximum_radius_correction < modification_tolerance
    && maximum_state_change < modification_tolerance
    && maximum_mass_error < invariant_tolerance
    && maximum_entropy_error < invariant_tolerance;

  std::cout << "zones=" << zones
            << " initial_residual=" << initial_residual
            << " final_residual=" << final_residual
            << " max_dr_over_width=" << maximum_radius_correction
            << " max_state_change=" << maximum_state_change
            << " max_shell_mass_error=" << maximum_mass_error
            << " max_entropy_error=" << maximum_entropy_error
            << " finite_ordered=" << finite_and_ordered
            << " passed=" << passed << '\n';
  return passed;
}

bool check_production_initializer(const bool yiming) {
  ThreeFluidSim sim;
  const int requested_zones = yiming ? 150 : 80;
  sim.initSolver(requested_zones);
  sim.initCoeffsYiming();
  if (yiming) {
    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);
  } else {
    sim.initPlummer(0.5, 1e-10, 1.0, 1.0, 1.0);
  }

  const auto initial_radius = sim.R;
  const auto initial_density = sim.Rho;
  const auto initial_pressure = sim.P;
  const auto initial_energy = sim.U;
  double initial_residual = 0.0;
  bool active_outer_shell = true;
  for (int f = 0; f < NF; ++f) {
    initial_residual = std::max(initial_residual,
                                hydrostatic_residual(sim, f));
    active_outer_shell = active_outer_shell
      && finite_bits(sim.Rho[f][sim.param.N-1])
      && finite_bits(sim.P[f][sim.param.N-1])
      && finite_bits(sim.U[f][sim.param.N-1])
      && sim.Rho[f][sim.param.N-1] > 0.0
      && sim.P[f][sim.param.N-1] > 0.0
      && sim.U[f][sim.param.N-1] > 0.0;
  }

  for (int f = 0; f < NF; ++f) {
    sim.solveRelaxationLAPACKE(f);
    sim.solveRelaxationLAPACKE(f);
  }
  sim.updateEnclosedMass();

  double relaxed_residual = 0.0;
  double maximum_relaxation_change = 0.0;
  for (int f = 0; f < NF; ++f) {
    relaxed_residual = std::max(relaxed_residual,
                                hydrostatic_residual(sim, f));
    for (int i = 0; i < sim.param.N; ++i) {
      maximum_relaxation_change = std::max({
        maximum_relaxation_change,
        relative_difference(sim.R[f][i], initial_radius[f][i]),
        relative_difference(sim.Rho[f][i], initial_density[f][i]),
        relative_difference(sim.P[f][i], initial_pressure[f][i]),
        relative_difference(sim.U[f][i], initial_energy[f][i])});
    }
  }

  const std::array<double, NF> mass_before_realign{
    sim.Menc[FS][sim.param.N-1], sim.Menc[FB][sim.param.N-1],
    sim.Menc[FD][sim.param.N-1]};
  sim.realign();
  double realigned_residual = 0.0;
  double maximum_realign_mass_error = 0.0;
  for (int f = 0; f < NF; ++f) {
    realigned_residual = std::max(realigned_residual,
                                  hydrostatic_residual(sim, f));
    maximum_realign_mass_error = std::max(
      maximum_realign_mass_error,
      relative_difference(sim.Menc[f][sim.param.N-1], mass_before_realign[f]));
    active_outer_shell = active_outer_shell
      && sim.Rho[f][sim.param.N-1] > 0.0 && sim.P[f][sim.param.N-1] > 0.0
      && sim.U[f][sim.param.N-1] > 0.0;
  }

  const std::array<double, NF> outer_energy_before{
    sim.U[FS][sim.param.N-1], sim.U[FB][sim.param.N-1], sim.U[FD][sim.param.N-1]};
  sim.solveConductionLAPACKE();
  double maximum_outer_conduction_change = 0.0;
  bool finite_after_conduction = true;
  for (int f = 0; f < NF; ++f) {
    maximum_outer_conduction_change = std::max(
      maximum_outer_conduction_change,
      relative_difference(sim.U[f][sim.param.N-1], outer_energy_before[f]));
    finite_after_conduction = finite_after_conduction
      && sim.U[f].array().isFinite().all()
      && sim.P[f].array().isFinite().all();
  }

  const bool passed = active_outer_shell && finite_after_conduction
    && initial_residual < residual_tolerance
    && relaxed_residual < residual_tolerance
    && realigned_residual < residual_tolerance
    && maximum_relaxation_change < modification_tolerance
    && maximum_realign_mass_error < invariant_tolerance
    && maximum_outer_conduction_change < invariant_tolerance;

  std::cout << "initializer=" << (yiming ? "Yiming" : "Plummer")
            << " zones=" << sim.param.N
            << " initial_residual=" << initial_residual
            << " relaxed_residual=" << relaxed_residual
            << " realigned_residual=" << realigned_residual
            << " max_relaxation_change=" << maximum_relaxation_change
            << " max_realign_mass_error=" << maximum_realign_mass_error
            << " max_outer_conduction_change="
            << maximum_outer_conduction_change
            << " active_outer_shell=" << active_outer_shell
            << " finite_after_conduction=" << finite_after_conduction
            << " passed=" << passed << '\n';
  return passed;
}

}  // namespace

int main() {
  std::cout << std::setprecision(17);
  std::cout << "residual_tolerance=" << residual_tolerance
            << " modification_tolerance=" << modification_tolerance
            << " invariant_tolerance=" << invariant_tolerance << '\n';
  bool passed = true;
  for (const int zones : {32, 80, 150}) {
    passed = check_grid(zones) && passed;
  }
  passed = check_production_initializer(false) && passed;
  passed = check_production_initializer(true) && passed;
  std::cout << "hydrostatic_relaxation_check_passed=" << passed << '\n';
  return passed ? 0 : 1;
}
