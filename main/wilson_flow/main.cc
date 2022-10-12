#include<Grid/Grid.h>
#include <common.h>

using namespace GridMeas;
using namespace Grid;

struct MeasArgs: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MeasArgs,
				  ActionType, action,
				  RealD, mass,
				  RealD, mobius_scale,
				  int, Ls,
				  std::vector<Integer>, GparityDirs,
				  std::string, cfg_stub,
				  std::string, rng_stub,
				  RealD, prop_invert_tol,
				  RealD, prop_invert_mixcg_innertol,
				  RealD, wflow_init_epsilon,
				  RealD, wflow_maxTau,
				  RealD, wflow_tolerance,
				  int, tslice_topq_meas_freq,
				  int, tslice_plaq_meas_freq);
  MeasArgs() { 
    action = ActionType::DWF;
    mass = 0.01;
    mobius_scale = 2.;
    Ls = 12;
    prop_invert_tol = 1e-8;
    prop_invert_mixcg_innertol = 1e-5;
    GparityDirs = {1,0,0};
    cfg_stub = "ckpoint_lat";
    rng_stub = "ckpoint_rng";
    wflow_init_epsilon = 0.01;
    wflow_maxTau = 8;
    wflow_tolerance = 1e-5;
    tslice_topq_meas_freq=20; //timeslice topq
    tslice_plaq_meas_freq=20; //timeslice plaq
  }
};

//Line format:
//tau Q[0] Q[1] ... Q[Lt-1]
//where tau is the Wilson flow time (or timeslice plaquette)
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
  
  Grids gridsD = makeDoublePrecGrids(args.Ls, latt);
  Grids gridsF = makeSinglePrecGrids(args.Ls, latt);

  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = args.GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams action_params;
  for(int i=0;i<3;i++) action_params.twists[i] = args.GparityDirs[i];
  action_params.twists[3] = 1; //antiperiodic BC in time

  LatticeGaugeFieldD U_d(gridsD.UGrid);
  LatticeGaugeFieldF U_f(gridsF.UGrid);
  LatticeGaugeFieldD V_d(gridsD.UGrid); //smeared

  Actions actions(args.action, action_params, args.mass, args.mobius_scale, U_d, gridsD, U_f, gridsF);

  //Start traj loop
  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    std::cout << GridLogMessage << "Starting traj " << traj << std::endl;

    std::vector<int> seeds4({traj, traj+2, traj+3, traj+4});
    GridParallelRNG pRNG(gridsD.UGrid); //4D!
    pRNG.SeedFixedIntegers(seeds4);
    
    GridSerialRNG sRNG;  
    sRNG.SeedFixedIntegers(seeds4); 

    std::cout << GridLogMessage << "Reading configuration" << std::endl;
    cps_cfg ? 
      readCPSconfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub) :
      readConfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);

    precisionChange(U_f,U_d);
    actions.ImportGauge(U_d,U_f);

    //Timeslice plaquette
    auto tslice_plaq = WilsonLoops<ConjugateGimplD>::timesliceAvgSpatialPlaquette(U_d);
    asciiWriteArray(tslice_plaq, "timeslice_plaq", traj);   
    
    //Chiral condensate
    {
      FermionFieldD rand_vol_src = randomGaussianVolumeSource(pRNG, gridsD.UGrid);
      FermionFieldD rand_vol_sol(gridsD.UGrid);
      mixedPrecInvertFieldXconj(rand_vol_sol, rand_vol_src, *actions.xconj_action_d, *actions.xconj_action_f, args.prop_invert_tol, args.prop_invert_mixcg_innertol);
      RealD cc = chiralCondensate(rand_vol_sol, rand_vol_src);
      asciiWriteValue(cc, "chiral_condensate", traj);
    }

    //Setup measurements to perform during the smearing
    WilsonFlowIO wflow_io;
    wflow_io.do_energy_density_clover = true;
    wflow_io.do_energy_density_plaq = true;
    wflow_io.do_timeslice_topq = true;
    wflow_io.do_timeslice_plaq = true;
    wflow_io.timeslice_topq_meas_freq = args.tslice_topq_meas_freq;
    wflow_io.timeslice_plaq_meas_freq = args.tslice_plaq_meas_freq;
    
    std::cout << GridLogMessage << "Starting Wilson Flow measurement" << std::endl;
    WilsonFlowAdaptiveMeasGeneral(wflow_io, args.wflow_init_epsilon, args.wflow_maxTau, args.wflow_tolerance, U_d, &V_d);
	
    asciiWriteArray(wflow_io.energy_density_clover, "wflow_clover", traj);
    asciiWriteArray(wflow_io.energy_density_plaq, "wflow_plaq", traj);
    writeTsliceTopQsmr(wflow_io.timeslice_topq, "timeslice_topq5li_smr", traj);
    writeTsliceTopQsmr(wflow_io.timeslice_plaq, "timeslice_plaq_smr", traj);

    //Measure topological charge
    std::vector<std::vector<Real> > tslice_topq5li_contribs = timesliceTopologicalCharge5LiContributions(V_d);
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
