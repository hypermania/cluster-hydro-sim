#include "three_fluid.hpp"
#include "observer.hpp"
#include "statler_reproduction.hpp"
#include "statler_observer.hpp"
#include "moving_comparison.hpp"
#include "moving_statler_observer.hpp"
#include "moving_statler_initialization.hpp"
#include "heggie_reproduction.hpp"
#include <filesystem>

void statler_reproduction(bool direct, const std::string& directory,
                          long long max_steps, double final_time_trh) {
  if(directory.empty()) throw std::invalid_argument("output directory is empty");
  if(!std::isfinite(final_time_trh) || final_time_trh<0)
    throw std::invalid_argument("invalid final time");
  const StatlerInitParam initial;
  const auto observing=statlerObserverParameters(initial);
  ThreeFluidSim sim;
  sim.param=statlerParameters();
  sim.param.maxSteps=max_steps;
  sim.param.maxTime=final_time_trh/observing.time_unit_over_trh;
  initializeCaptureCluster(sim,initial);
  if(!direct) sim.param.c4.fill(0);
  std::filesystem::create_directories(directory);
  const std::string output=directory.back()=='/'?directory:directory+"/";
  sim.saveParams(output);
  std::filesystem::create_directories(output+"initialization");
  std::filesystem::create_directories(output+"observer");
  save_param_for_Mathematica(initial,output+"initialization/");
  save_param_for_Mathematica(observing,output+"observer/");
  StatlerObserver observer(observing);
  sim.evolve(observer);
  observer.save(output);
}

void moving_statler_reproduction(const std::string& directory,int zones,
                                double final_time_trh,long long max_steps) {
  if(directory.empty()||zones<3||!std::isfinite(final_time_trh)||final_time_trh<=0||max_steps<=0)
    throw std::invalid_argument("invalid moving Statler controls");
  if(std::filesystem::exists(directory+"/param.dat"))
    throw std::invalid_argument("choose a fresh output directory");
  const StatlerInitParam initial;
  auto observing=statlerObserverParameters(initial);
  // Approximate square epochs read from Statler Fig. 11a; retain the actual
  // accepted times in every legend. The paper does not tabulate these times.
  observing.snapshot_times_trh={0,1,5,10,18,50,350,2000,10000};
  MovingStatlerGrid grid;grid.zones=zones;
  MovingThreeFluidSim sim;
  sim.param.Deltat=1e-9;sim.param.max_timestep=1e4;
  sim.param.maxTime=final_time_trh/observing.time_unit_over_trh;
  sim.param.maxSteps=max_steps;
  initializeMovingCaptureCluster(sim,initial,grid);
  const std::string output=directory.back()=='/'?directory:directory+"/";
  std::filesystem::create_directories(output+"initialization");
  std::filesystem::create_directories(output+"observer");
  std::filesystem::create_directories(output+"grid");
  save_param_for_Mathematica(initial,output+"initialization/");
  save_param_for_Mathematica(observing,output+"observer/");
  save_param_for_Mathematica(grid,output+"grid/");
  save_param_for_Mathematica(sim.param,output);
  write_to_file(sim.faces(),output+"initial_faces.dat");
  MovingStatlerObserver observer(observing,output);
  try {sim.evolve(observer);}catch(...) {observer.save(output);throw;}
  observer.save(output);
  if(sim.totalTime<sim.param.maxTime)
    throw std::runtime_error("moving Statler run stopped before requested endpoint");
}

void one_fluid_split_in_two(void){
  const std::vector times_to_save({0.0, 5.5, 5.565, 5.56517});
  // const std::vector times_to_save({0.0, 0.1, 1.0, 5.5, 5.53, 5.538});

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.param.md = sim.param.ms;
    sim.param.c2[FD] = sim.param.c2[FS];

    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/one_fluid_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.param.md = sim.param.ms;
    sim.param.c2[FD] = sim.param.c2[FS];
  
    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/one_fluid_split_in_two_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
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
    sim.param.md = sim.param.ms;
    sim.param.c2[FD] = sim.param.c2[FS];

    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);
    
    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);
    
    std::string dir = "output/one_fluid_without_tidal/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.param.md = sim.param.ms;
    sim.param.c2[FD] = sim.param.c2[FS];
    sim.param.tidal_cutoff = 1;
    sim.param.tidal_radius = 10.0;
    sim.param.tidal_cutoff_factor = 50.0;

    // sim.printParams();
    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);

    std::string dir = "output/one_fluid_with_tidal/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
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
    sim.param.md = 0.1 * sim.param.ms;
    // sim.param.md = sim.param.ms;
    // sim.param.c2[FD] = 0.1 * sim.param.c2[FS];

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    //LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/baseline_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually set DM mass
    sim.param.md = 0.1 * sim.param.ms;
    sim.param.tidal_cutoff = 1;
    sim.param.tidal_cutoff_factor = 10;
    sim.param.tidal_radius = 2;

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    // LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ApproximateTimeObserver observer1({0.0, 1.0, 1.8, 1.89, 1.898});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/with_tidal_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
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
    const double Mtot = sim.Menc[FS][sim.param.N-1] + sim.Menc[FB][sim.param.N-1] + sim.Menc[FD][sim.param.N-1];
    sim.param.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.param.mb = 2.0 * sim.param.ms;
    sim.param.md = 1e-10 * sim.param.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    // // Artificially turn off dynamical heating to check
    // for(int f = 0; f < NF; ++f){
    //   for(int f2 = 0; f2 < NF; ++f2){
    // 	sim.param.c1[f*NF+f2] = 0.0;
    //   }
    // }

    ApproximateTimeObserver observer1({0.0, 5.0, 5.8, 5.806, 5.80664});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);

    
    std::string dir = "output/one_fluid_baseline/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);
  }
  
  // With binary formation, no heating
  {
    ThreeFluidSim sim;
    sim.initSolver(N);
    sim.initPlummer(1.0, 1e-10, 1e-10, 1.0, 1.0);
    const double Mtot = sim.Menc[FS][sim.param.N-1] + sim.Menc[FB][sim.param.N-1] + sim.Menc[FD][sim.param.N-1];
    sim.param.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.param.mb = 2.0 * sim.param.ms;
    sim.param.md = 1e-10 * sim.param.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    sim.param.binary_formation = BINARY_FORMATION_MODE_2;
    // Artificially turn off binary heating to check
    for(int f = 0; f < NF; ++f){
      for(int f2 = 0; f2 < NF; ++f2){
	sim.param.c4[(f)*NF+(f2)] = 0.0;
      }
    }

    ApproximateTimeObserver observer1({0.0, 3.0, 3.8, 3.82, 3.823});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);    

    
    std::string dir = "output/one_fluid_binary_formation_without_binary_heating/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);

    std::cout << "Menc[FS], Menc[FB] = "
	      << sim.Menc[FS][N-1] << "," << sim.Menc[FB][N-1] << std::endl;
  }
  
  // With binary formation and heating
  {
    ThreeFluidSim sim;
    sim.initSolver(N);
    sim.initPlummer(1.0, 1e-10, 1e-10, 1.0, 1.0);
    const double Mtot = sim.Menc[FS][sim.param.N-1] + sim.Menc[FB][sim.param.N-1] + sim.Menc[FD][sim.param.N-1];
    sim.param.ms = Mtot / 1e6; // Assuming cluster mass = 10^6 ms
    sim.param.mb = 2.0 * sim.param.ms;
    sim.param.md = 1e-10 * sim.param.ms;
    sim.initCoeffs(1e6, 2.0, 1e-10);
    sim.param.binary_formation = BINARY_FORMATION_MODE_2;

    ApproximateTimeObserver observer1({0.0, 1.0, 10.0, 100.0, 692.0});
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    KeyValueObserver observer3;
    ObserverPack observer(observer1, observer2, observer3);    
    
    std::string dir = "output/one_fluid_binary_formation/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(observer);
    observer.save(dir);

    std::cout << "Menc[FS], Menc[FB] = "
	      << sim.Menc[FS][N-1] << "," << sim.Menc[FB][N-1] << std::endl;
  }

}


int main(int argc, char** argv) {
  if(argc>1) {
    try {
      if(std::string(argv[1])=="moving-heggie") {
        if(argc<4||argc>6)throw std::invalid_argument(
          "usage: main moving-heggie output_directory model(0:single,1:segregation,2:BS,3:BS+BB) [f_star] [zones]");
        HeggieInitParam initial;initial.model=std::stoll(argv[3]);
        if(argc>4)initial.f_star=std::stod(argv[4]);
        if(argc>5)initial.zones=std::stoll(argv[5]);
        runMovingHeggie(argv[2],initial);return 0;
      }
      if(std::string(argv[1])=="moving-statler") {
        if(argc<3||argc>6)throw std::invalid_argument(
          "usage: main moving-statler output_directory [zones=500] [final_time_trh=10000] [max_steps=2000000]");
        moving_statler_reproduction(argv[2],argc>3?std::stoi(argv[3]):500,
          argc>4?std::stod(argv[4]):10000,argc>5?std::stoll(argv[5]):2000000);
        return 0;
      }
      if(std::string(argv[1])=="moving-comparison") {
        if(argc<5||argc>9)throw std::invalid_argument(
          "usage: main moving-comparison output_directory sample(0..3) final_time [epsilon] [max_dt] [canonical_zones] [max_steps]");
        runMovingComparison(argv[2],std::stoi(argv[3]),std::stod(argv[4]),
                            argc>5?std::stod(argv[5]):0,argc>6?std::stod(argv[6]):1e-3,
                            argc>7?std::stoi(argv[7]):500,argc>8?std::stoll(argv[8]):20000);
        return 0;
      }
      if(argc<4 || argc>6 || std::string(argv[1])!="statler" ||
         (std::string(argv[2])!="direct" && std::string(argv[2])!="control"))
        throw std::invalid_argument("usage: main statler direct|control output_directory [max_steps] [final_time_trh]");
      statler_reproduction(std::string(argv[2])=="direct",argv[3],
                            argc>4?std::stoll(argv[4]):2000000,
                            argc>5?std::stod(argv[5]):10000.0);
      return 0;
    } catch(const std::exception& e) {
      std::cerr<<e.what()<<'\n';return 1;
    }
  }
  // one_fluid_split_in_two();
  // tidal_bench_single();
  // tidal_bench_AB();
  binary_formation();
  
  return 0;
}
