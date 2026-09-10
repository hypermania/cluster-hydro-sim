#include "three_fluid.hpp"
#include "observer.hpp"

void one_fluid_split_in_two(void){
  const std::vector times_to_save({0.0, 5.5, 5.565, 5.56517});
  // const std::vector times_to_save({0.0, 0.1, 1.0, 5.5, 5.53, 5.538});

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];

    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/one_fluid_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];
  
    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/one_fluid_split_in_two_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

}

void tidal_bench_single(void){
  const std::vector times_to_save({0.0, 0.1, 1.0, 5.5, 5.53, 5.538});

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];

    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);
    
    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);
    
    std::string dir = "output/one_fluid_without_tidal/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];
    sim.tidal_cutoff = 1;
    sim.tidal_radius = 10.0;
    sim.tidal_cutoff_factor = 50.0;

    // sim.printParams();
    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);

    std::string dir = "output/one_fluid_with_tidal/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(100000, observer);
    observer.save(dir);
  }
  
}

void tidal_bench_AB(void){
  const std::vector times_to_save({0.0, 3.0, 3.2, 3.24, 3.243});
  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually set DM mass
    sim.md = 0.1 * sim.ms;
    // sim.md = sim.ms;
    // sim.c2[FD] = 0.1 * sim.c2[FS];

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    //LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/baseline_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually set DM mass
    sim.md = 0.1 * sim.ms;
    sim.tidal_cutoff = 1;
    sim.tidal_cutoff_factor = 10;
    sim.tidal_radius = 2;

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    // LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ApproximateTimeObserver observer1({0.0, 1.0, 1.8, 1.89, 1.898});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/with_tidal_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

}

void binary_formation(void){
  const int N = 500;
  
  // Baseline
  {
    ThreeFluidSim sim;
    sim.initSolver(N);
    sim.initPlummer(1.0, 1e-10, 1e-10, 1.0, 1.0);
    const double Mtot = sim.Menc[FS][sim.N-1] + sim.Menc[FB][sim.N-1] + sim.Menc[FD][sim.N-1];
    sim.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.mb = 2.0 * sim.ms;
    sim.md = 1e-10 * sim.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    // // Artificially turn off dynamical heating to check
    // for(int f = 0; f < NF; ++f){
    //   for(int f2 = 0; f2 < NF; ++f2){
    // 	sim.c1[f][f2] = 0.0;
    //   }
    // }

    ApproximateTimeObserver observer1({0.0, 5.0, 5.8, 5.806, 5.80664});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);

    
    std::string dir = "output/one_fluid_baseline/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }
  
  // With binary formation, no heating
  {
    ThreeFluidSim sim;
    sim.initSolver(N);
    sim.initPlummer(1.0, 1e-10, 1e-10, 1.0, 1.0);
    const double Mtot = sim.Menc[FS][sim.N-1] + sim.Menc[FB][sim.N-1] + sim.Menc[FD][sim.N-1];
    sim.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.mb = 2.0 * sim.ms;
    sim.md = 1e-10 * sim.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    sim.binary_formation = BINARY_FORMATION_MODE_2;
    // Artificially turn off binary heating to check
    for(int f = 0; f < NF; ++f){
      for(int f2 = 0; f2 < NF; ++f2){
	sim.c4[f][f2] = 0.0;
      }
    }

    ApproximateTimeObserver observer1({0.0, 3.0, 3.8, 3.82, 3.823});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);    

    
    std::string dir = "output/one_fluid_binary_formation_without_binary_heating/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(1000000, observer);
    observer.save(dir);

    std::cout << "Menc[FS], Menc[FB] = "
	      << sim.Menc[FS][N-1] << "," << sim.Menc[FB][N-1] << std::endl;
  }
  
  // With binary formation and heating
  {
    ThreeFluidSim sim;
    sim.initSolver(N);
    sim.initPlummer(1.0, 1e-10, 1e-10, 1.0, 1.0);
    const double Mtot = sim.Menc[FS][sim.N-1] + sim.Menc[FB][sim.N-1] + sim.Menc[FD][sim.N-1];
    sim.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.mb = 2.0 * sim.ms;
    sim.md = 1e-10 * sim.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    sim.binary_formation = BINARY_FORMATION_MODE_2;

    ApproximateTimeObserver observer1({0.0, 1.0, 10.0, 100.0, 692.0});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);    
    
    std::string dir = "output/one_fluid_binary_formation/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(1000000, observer);
    observer.save(dir);

    std::cout << "Menc[FS], Menc[FB] = "
	      << sim.Menc[FS][N-1] << "," << sim.Menc[FB][N-1] << std::endl;
  }

}


int main() {
  // one_fluid_split_in_two();
  // tidal_bench_single();
  // tidal_bench_AB();
  binary_formation();
  
  return 0;
}
