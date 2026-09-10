#pragma once

#include "Eigen/Dense"
#include <vector>
#include "param.hpp"
#include "io.hpp"



constexpr double PI = 3.14159265358979323846;
constexpr double C1 = 1.0 / sqrt(3.0 * PI);          // ≈0.326
constexpr double ALPHA = 2 * erf(sqrt(1.5)) - 2 * exp(-1.5) * sqrt(6.0 / PI);
constexpr double BETA_C = 0.45;
constexpr double C2 = sqrt(2.0) * pow(3.0, -2) * ALPHA * BETA_C; // ≈0.086
constexpr double GAMMA_LAMBDA = 0.4;
constexpr double GAMMA_COULOMB = 0.11;

// Fluid indices (avoid name clash with any variable)
constexpr int FS = 0;  // single stars
constexpr int FB = 1;  // binaries
constexpr int FD = 2;  // dark matter
constexpr int NF = 3;

// Simulation flags
constexpr long long int BINARY_FORMATION_OFF = 0;
constexpr long long int BINARY_FORMATION_MODE_1 = 1;
constexpr long long int BINARY_FORMATION_MODE_2 = 2;
constexpr long long int BINARY_FORMATION_POWER_LAW = 3;

constexpr long long int TIDAL_CUTOFF_OFF = 0;
constexpr long long int TIDAL_CUTOFF_ON = 1;


struct ThreeFluidParam {
  // Evolution parameters
  long long int N;
  double ms;
  double mb;
  double md;
  std::array<double, NF> c2; // Conduction
  std::array<double, NF*NF> c1; // Dynamical heating
  std::array<double, NF*NF> c4;  // Binary heating
  
  long long int binary_formation;
  long long int tidal_cutoff;
  double tidal_cutoff_factor;
  double tidal_radius;

  // Numerical control parameters
  double Deltat;
  double StopDensity;
  double maxTime;
  long long int maxSteps;
  double thres;
  // Additional configuration. Physical scales are metadata, not internal state.
  long long int schema_version = 2;
  long long int initial_profile = 0; // 0 Plummer, 1 Yiming
  std::array<double, 5> initial_profile_values = {1, 1e-10, 1e-10, 1, 1};
  double initial_timestep = 1e-3;
  long long int timestep_controller = 0; // 0 legacy U, 1 bounded U+density
  double max_timestep = 1.0;
  double donor_fraction_limit = 0.005;
  double retry_safety = 0.8;
  long long int max_retries = 100;
  double density_change_tolerance = 0.025;
  double timestep_growth_min = 0.5;
  double timestep_growth_max = 1.5;
  double error_control_floor = 1e-12;
  long long int relaxation_passes = 2;
  long long int central_density_measure = 0; // 0 max species, 1 singles+binaries
  double capture_coefficient = 0; // number source: A rho_s^2 U_s^-exponent
  double capture_energy_exponent = 0.6;
  double captured_specific_energy_fraction = 0.5;
  double reference_coulomb_log = 0;
  double stellar_number = 0;
  double stellar_mass_msun = 0;
  double stellar_radius_rsun = 0;
  double length_unit_pc = 0;
  double reference_trh_myr = 0;
  double gravitational_constant_pc3_msun_myr2 = 0;
  double solar_radius_pc = 0;
  double time_unit_myr = 0;
  double time_unit_over_trh = 0;
  double final_time_trh = 0;
  long long int direct_heating = 1;
  // Observer settings: explicit finite schedule; unused entries are zero.
  long long int statler_observer = 0;
  long long int history_stride = 1;
  long long int snapshot_count = 0;
  std::array<double, 16> snapshot_times_trh{};
  long long int save_peak_snapshot = 1;
  long long int output_format = 1; // little-endian binary float64/int64
};


// ===================================================================
// Full 3-fluid simulation
// All units are dimensionless, in terms of arbitrary r0, rho0, t0,
// Mass unit is in M0 = 4 pi rho0 r0^3
// See arXiv:2505:18251, eqn (35)-(39)
// ===================================================================
class ThreeFluidSim {
public:
  // Evolution parameters
  long long int N = 200;
  double ms = 1.0 / 1e6;
  double mb = 2.0 / 1e6;
  double md = 1e-10 / 1e6;
  double c2[NF]; // Conduction
  double c1[NF][NF]; // Dynamical heating
  double c4[NF][NF]; // Binary heating

  long long int binary_formation = BINARY_FORMATION_OFF;
  long long int tidal_cutoff = TIDAL_CUTOFF_OFF;
  double tidal_cutoff_factor = 50;
  double tidal_radius = 1e1;


  // Numerical control parameters
  double Deltat = 1e-3;
  double StopDensity = 1e12;
  double maxTime = 1e4;
  long long int maxSteps = 1e7;
  double thres = 1e-3;
  
  // New settings live here; parameters() exports these plus legacy members.
  ThreeFluidParam options{};
  double cumulative_formed_binaries = 0.0;
  long long int rejected_steps = 0;
  std::array<Eigen::VectorXd, NF> trialU, trialP;
  Eigen::VectorXd formationSource;

  // Internal state
  double totalTime = 0.0;
  long long int step = 0;
  std::array<Eigen::VectorXd, NF> Rho, U, Menc, P, R;

  // Temporary variables
  std::vector<double> conductionAB;
  std::vector<double> conductionBSTD;
  std::vector<int> conductionIPIV;
  
  std::array<Eigen::VectorXd, NF> sqrtU;
  std::array<Eigen::VectorXd, NF> logRho;
  Eigen::VectorXd logR;

  Eigen::VectorXd hydroDL;
  Eigen::VectorXd hydroD;
  Eigen::VectorXd hydroDU;
  Eigen::VectorXd hydroB;

  Eigen::VectorXd newR;
  Eigen::VectorXd newRho;

  // Initialize empty solver
  ThreeFluidSim();

  // Allocate and zero memory for the solver
  void initSolver(const int N);

  // Reset time evolution
  void initControl();
  
  // Compute and assign coeffs c1[NF][NF], c2[NF] and c4[NF][NF]
  void initCoeffs(const double Mtot_over_ms, const double mb_over_ms, const double md_over_ms);

  // Assign coeffs according to Yiming's paper
  void initCoeffsYiming();

  // Assign initial conditions using algorithm replicated from Yiming's paper
  void initPlummerYiming(const double rho0, const double xi1, const double xi2, const double zeta1, const double zeta2);
  // Assign initial conditions
  void initPlummer(const double rhos_central, const double xi1, const double xi2, const double zeta1, const double zeta2);

  // IO
  void printParams() const;
  void printCoeffs() const;
  void saveParams(const std::string &) const;
  ThreeFluidParam parameters() const;


  // Fragments of main evolution step
  void updateEnclosedMass();
  void solveConductionLAPACKE();
  void solveRelaxationLAPACKE(const int f);
  void realign();
  void applyBinaryFormation();
  void applyTidalCutoff();
  double centralDensity() const;
  double captureNumberRate(int zone) const;
  double prepareFormationSource();
  void applyPowerLawFormation();
  void projectHydrostatic();
  void validateEvolution() const;
  void advanceAcceptedStep();
  double conductionChange() const;
  void selectNextTimestep(double used_dt, double energy_change, double old_density);
  bool stopCondition() const;
  void sanityCheck() const;

  // Main evolution step
  template<typename Observer>
  void evolve(Observer &observer) {
    validateEvolution();
    step = 0;
    observer(static_cast<const ThreeFluidSim&>(*this));
    while (!stopCondition()) {
      advanceAcceptedStep();
      observer(static_cast<const ThreeFluidSim&>(*this));
    }
  }
};

