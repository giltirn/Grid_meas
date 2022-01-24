#include<Grid/Grid.h>
#include <common.h>

using namespace GridMeas;
using namespace Grid;

struct MeasArgs: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MeasArgs,
				  std::vector<Integer>, GparityDirs,
				  std::string, cfg_stub,
				  std::string, rng_stub,
				  int, Nstep,
				  RealD, epsilon,
				  int, topq_meas_freq);
  MeasArgs() { 
    GparityDirs = {1,0,0};
    cfg_stub = "ckpoint_lat";
    rng_stub = "ckpoint_rng";
    Nstep = 100;
    epsilon = 0.01;
    topq_meas_freq=20;
  }
};

//Line format:
//tau Q[0] Q[1] ... Q[Lt-1]
//where tau is the Wilson flow time 
void writeTsliceTopQsmr(const std::vector<std::pair<RealD,std::vector<RealD> > > &data,
			const std::string &stub, const int traj){
    std::string filename = stub + "." + std::to_string(traj);
    std::ofstream of(filename);
    of.precision(17);

    for(int t=0;t<data.size();t++){
      of << data[t].first;
      for(int tt=0;tt<data[t].second.size();tt++)
	of << " " << data[t].second[tt];
      of << std::endl;
    }
}

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  assert(argc >= 5);
  std::string arg_file = argv[1];
  int cfg_start = std::stoi(argv[2]);
  int cfg_lessthan = std::stoi(argv[3]);
  int cfg_step = std::stoi(argv[4]);

  bool cps_cfg = false;
  for(int i=1;i<argc;i++){
    std::string sargv(argv[i]);
    if(sargv == "--cps_cfg"){
      cps_cfg = true;
    }
  }

  MeasArgs args;
  
  if(fileExists(arg_file)){
    std::cout << GridLogMessage << " Reading " << arg_file << std::endl;
    Grid::XmlReader rd(arg_file);
    read(rd, "Args", args);
  }else if(!GlobalSharedMemory::WorldRank){
    std::cout << GridLogMessage << " File " << arg_file << " does not exist" << std::endl;
    std::cout << GridLogMessage << " Writing xml template to " << arg_file << ".templ" << std::endl;
    Grid::XmlWriter wr(arg_file + ".templ");
    write(wr, "Args", args);

    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }

  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];

  auto UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  auto UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = args.GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  LatticeGaugeFieldD U(UGrid);
  LatticeGaugeFieldD V(UGrid); //smeared

  //Start traj loop
  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    std::cout << GridLogMessage << "Starting traj " << traj << std::endl;

    std::vector<int> seeds4({traj, traj+2, traj+3, traj+4});
    GridParallelRNG pRNG(UGrid); //4D!
    pRNG.SeedFixedIntegers(seeds4);
    
    GridSerialRNG sRNG;  
    sRNG.SeedFixedIntegers(seeds4); 

    std::cout << GridLogMessage << "Reading configuration" << std::endl;
    cps_cfg ? 
      readCPSconfiguration(U, sRNG, pRNG, traj, args.cfg_stub) :
      readConfiguration(U, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);

    //Timeslice plaquette
    auto tslice_plaq = WilsonLoops<ConjugateGimplD>::timesliceAvgSpatialPlaquette(U);
    asciiWriteArray(tslice_plaq, "timeslice_plaq", traj);   

    //Do Wilson flow and measure topological charge along the way
    std::vector<std::pair<RealD,RealD> > wflow;
    std::vector<std::pair<RealD,std::vector<RealD> > > tslice_topq5li_int; //Q[t] during smearing

    std::cout << GridLogMessage << "Starting Wilson Flow measurement" << std::endl;
    WilsonFlowEnergyDensityAndTimesliceTopQ(wflow, tslice_topq5li_int, args.Nstep, args.epsilon, args.topq_meas_freq, U, &V);

    //auto wflow = WilsonFlowEnergyDensity(args.Nstep, args.epsilon, U, &V);
    asciiWriteArray(wflow, "wflow", traj);
    writeTsliceTopQsmr(tslice_topq5li_int, "timeslice_topq5li_smr", traj);

    //Measure topological charge
    std::vector<std::vector<Real> > tslice_topq5li_contribs = timesliceTopologicalCharge5LiContributions(V);
    asciiWriteArray(tslice_topq5li_contribs, "timeslice_topq5li_contribs", traj);

    std::vector<Real> tslice_topq5li = timesliceTopologicalCharge5Li(tslice_topq5li_contribs);
    asciiWriteArray(tslice_topq5li, "timeslice_topq5li", traj);

    std::vector<Real> topq5li_contribs = topologicalCharge5LiContributions(tslice_topq5li_contribs);
    asciiWriteArray(topq5li_contribs, "topq5li_contribs", traj);
   
    RealD topq5li = topologicalCharge5Li(topq5li_contribs);
    asciiWriteValue(topq5li, "topq5li", traj);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
