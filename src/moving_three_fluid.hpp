#pragma once
#include "three_fluid.hpp"

constexpr long long HEATING_DONOR_DISPERSION=0, HEATING_RELATIVE_DISPERSION=1;

// Dimensionless evolution controls only. Initial profiles/grid and observers
// are supplied separately. Explicit stripping is not supported.
struct MovingThreeFluidParam {
  double epsilon = 1e-10;
  double q = 0;
  std::array<double, NF> mass{1e-6, 2e-6, 1e-16};
  std::array<double, NF> c2{};
  std::array<double, NF*NF> c1{}, c4{};
  long long binary_formation = BINARY_FORMATION_OFF;
  double capture_coefficient = 0; // number source A rho_s^2 U_s^-0.6
  long long heating_dispersion = HEATING_DONOR_DISPERSION;
  long long reflecting_boundary = 0; // insulated rigid wall; otherwise vacuum outflow
  double thermal_length_over_radius = 1;
  double Deltat = 1e-5;
  double max_timestep = 1e-3;
  double change_tolerance = 1e-3;
  double maxTime = 5;
  double StopDensity = 1e12;
  long long maxSteps = 1000000;
};

class MovingThreeFluidSim {
public:
  using Vec = Eigen::Vector3d;
  using Mat = Eigen::Matrix3d;
  static constexpr int KL=13, KU=14, LDAB=2*KL+KU+1;
  enum Field { RHO=0, VEL=1, U=2, MASS=3 };
  static constexpr int slot(int f, int field) {
    return field==MASS ? 9+f : 3*f+2-field;
  }
  static constexpr int index(int i,int f,int field) {return 12*i+slot(f,field);}
  MovingThreeFluidParam param;
  // One cell-major primitive state, followed by no hidden fluid arrays.
  std::vector<double> state;
  double totalTime=0, Deltat=0;
  long long step=0;
  std::array<double,NF> escaped_mass{}, escaped_energy{}, escaped_heat{};
  std::array<double,NF> last_mass_outflow{}, last_energy_outflow{}, last_heat_outflow{};
  double linear_residual=0;
  double last_change=0;
  int limiting_index=0;
  double last_gravity_work=0, last_heating=0, energy_ledger_error=0;
  double cumulative_formed_binaries=0, cumulative_capture_energy=0;
  double last_formed_binaries=0, last_capture_energy=0;

  // Input: positive profiles with a common mesh. reference_balance=true is
  // reserved for a supplied q=0 hydrostatic profile; q is never balanced away.
  void initialize(const std::vector<double>& faces,
                  const std::vector<double>& initial, bool reference_balance);
  int zones() const {return static_cast<int>(state.size()/12);}
  double value(int i,int f,int field) const {return state.at(index(i,f,field));}
  const std::vector<double>& faces() const {return edges;}
  double radius(int i) const {return centres.at(i);}
  double volume(int i) const {return volumes.at(i);}
  double centralDensity() const;
  void validate() const;
  void assembleStep();
  void advanceAcceptedStep();
  bool stopCondition() const;
  // Read-only assembled equations for independent dense/Jacobian tests.
  const std::vector<double>& bandMatrix() const {return band;}
  const std::vector<double>& rightHandSide() const {return rhs;}
  const std::vector<double>& lastIncrement() const {return increment;}
  template<class Observer> void evolve(Observer& observer) {
    if(step==0) Deltat=param.Deltat;
    validate();
    observer(static_cast<const MovingThreeFluidSim&>(*this));
    while(!stopCondition()) {
      advanceAcceptedStep();
      observer(static_cast<const MovingThreeFluidSim&>(*this));
    }
  }

private:
  struct Face {
    Vec flux=Vec::Zero();
    Mat left=Mat::Zero(), right=Mat::Zero();
    double heat=0, heat_left=0, heat_right=0;
  };
  std::vector<double> edges, centres, volumes, eta;
  std::vector<double> correction; // Frozen reference acceleration, not force density.
  std::vector<Face> fluxes;
  std::vector<double> band, rhs, factor, increment, scale, next, refinement;
  std::vector<int> pivots;
  Vec primitive(int i,int f) const;
  Vec storage(const Vec& x) const;
  Vec flux(const Vec& x) const;
  Mat storageJacobian(const Vec& x) const;
  Mat fluxJacobian(const Vec& x) const;
  Face boundaryFlux(const Vec& x) const;
  void thermalSource(int cell,Vec& source,Mat& derivative) const;
  void buildFluxes();
  void add(int row,int col,double value);
  void assembleCells();
  // Frozen per-donor capture frequency; local conservative transfer in the
  // same band solve. New binaries inherit v_s and random specific energy U_s/2.
  double captureFrequency(int cell) const;
  void formationSource(int cell,int fluid,Vec& source,Mat& derivative) const;
  void solveBanded();
  double recoverAndBudget();
  void selectNextTimestep(double change);
};
