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
  long long int N = 200;
  double ms = 1.0 / 1e6;
  double mb = 2.0 / 1e6;
  double md = 1e-10 / 1e6;
  std::array<double, NF> c2{}; // Conduction
  std::array<double, NF*NF> c1{}; // Dynamical heating, row-major
  std::array<double, NF*NF> c4{}; // Binary heating, row-major
  
  long long int binary_formation = BINARY_FORMATION_OFF;
  long long int tidal_cutoff = TIDAL_CUTOFF_OFF;
  double tidal_cutoff_factor = 50;
  double tidal_radius = 10;

  // Numerical control parameters
  double Deltat = 1e-3; // Initial timestep; adaptive timestep is simulation state.
  double StopDensity = 1e12;
  double maxTime = 1e4;
  long long int maxSteps = 10000000;
  double u_change_tolerance = 1e-3;
  double max_timestep = 1.0;
  double capture_coefficient = 0; // PT number source: A rho_s^2 U_s^-0.6
  long long int runtime_validation = 1; // Optional full-state scans.
};


// ===================================================================
// Full 3-fluid simulation
// All units are dimensionless, in terms of arbitrary r0, rho0, t0,
// Mass unit is in M0 = 4 pi rho0 r0^3
// See arXiv:2505:18251, eqn (35)-(39)
// ===================================================================
class ThreeFluidSim {
public:
  ThreeFluidParam param{};
  // Adaptive timestep is evolving state, not a second configuration copy.
  double Deltat = 1e-3;
  double cumulative_formed_binaries = 0.0;
  std::array<Eigen::VectorXd, NF> previousU; // Pre-conduction values for timestep control.

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

  // Deprecated compatibility initializer; initCoeffs is canonical.
  void initCoeffsYiming();

  // Deprecated compatibility initializer; initPlummer is canonical.
  void initPlummerYiming(const double rho0, const double xi1, const double xi2, const double zeta1, const double zeta2);
  // Assign initial conditions
  void initPlummer(const double rhos_central, const double xi1, const double xi2, const double zeta1, const double zeta2);

  // IO
  void printParams() const;
  void printCoeffs() const;
  void saveParams(const std::string &) const;


  // Fragments of main evolution step
  void updateEnclosedMass();
  void solveConductionLAPACKE();
  void solveRelaxationLAPACKE(const int f);
  void realign();
  void applyBinaryFormation();
  void applyTidalCutoff();
  double centralDensity() const;
  double captureNumberRate(int zone) const;
  void applyPowerLawFormation();
  void projectHydrostatic();
  void validateEvolution() const;
  void advanceAcceptedStep();
  double conductionChange() const;
  void selectNextTimestep(double used_dt, double energy_change);
  bool stopCondition() const;
  void sanityCheck() const;

  // Main evolution step
  template<typename Observer>
  void evolve(Observer &observer) {
    if(totalTime==0) Deltat=param.Deltat;
    validateEvolution();
    step = 0;
    observer(static_cast<const ThreeFluidSim&>(*this));
    while (!stopCondition()) {
      advanceAcceptedStep();
      observer(static_cast<const ThreeFluidSim&>(*this));
    }
  }
};
