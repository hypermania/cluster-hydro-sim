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


  // Fragments of main evolution step
  void updateEnclosedMass();
  void solveConductionLAPACKE();
  void solveRelaxationLAPACKE(const int f);
  void realign();
  void applyBinaryFormation();
  void applyTidalCutoff();
  bool stopCondition() const;
  void sanityCheck() const;

  // Main evolution step
  template<typename Observer>
  void evolve(Observer &observer) {
    using namespace std;
    using namespace Eigen;
    std::array<Eigen::VectorXd, NF> lastU(U);

    step = 0;
    while(true) {
      // std::cout << std::setprecision(9) << std::left;
      // cout << "step, t, Deltat = " << step << ", " << totalTime << ", " << Deltat << endl;
      observer(*this);
    
      if(stopCondition()) break;
      
      // Start of timestep
      lastU = U;

      solveConductionLAPACKE();

      // Adaptive timestep, should be calculated due to change from conduction step only
      double maxChange = 0.0;
      for(int f = 0; f < NF; ++f) {
	maxChange = max(maxChange, ((U[f].array() - lastU[f].array()).abs() / lastU[f].array()).maxCoeff());
      }

      if(tidal_cutoff != TIDAL_CUTOFF_OFF) { applyTidalCutoff(); }
      // Binary formation changes Rho and Menc, but preserves U, updates P = (2/3) * Rho * U
      if(binary_formation != BINARY_FORMATION_OFF) { applyBinaryFormation(); }

      
      for(int f = 0; f < NF; ++f) {
	solveRelaxationLAPACKE(f);
	solveRelaxationLAPACKE(f);
	// solveRelaxationLAPACKE(f);
	// solveRelaxationLAPACKE(f);
	// solveRelaxationLAPACKE(f);
      }
      
      realign();

      // Sanity checks
      sanityCheck();
      
      totalTime += Deltat;      
      Deltat = Deltat * thres / maxChange;
      ++step;
    }
    
    cout << "Simulation finished at t = " << totalTime << " (steps = " << step << ")" << endl;
  }

  
};


