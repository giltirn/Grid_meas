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
				  MixedCGargs, cg_args,
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

struct GparitySetup{
  typedef Actions ActionsType;
  typedef ConjugateGimplD GimplD;
  typedef CayleyFermion5D<GparityWilsonImplD> ActionTypeD;
  
  inline static std::unique_ptr<ActionsType> setupActions(const MeasArgs &args, LatticeGaugeFieldD &U_d, Grids &gridsD, LatticeGaugeFieldF &U_f, Grids &gridsF){ 
    GparityWilsonImplD::ImplParams action_params = setupActionParams(args.GparityDirs);
    return std::unique_ptr<Actions>(new Actions(args.action, action_params, args.mass, args.mobius_scale, U_d, gridsD, U_f, gridsF));
  }

  inline static void invertRandomSource(FermionFieldD &rand_vol_sol, const FermionFieldD &rand_vol_src, Actions &actions, const MixedCGargs &cg_args){
    mixedPrecInvertFieldXconj(rand_vol_sol, rand_vol_src, *actions.xconj_action_d, *actions.xconj_action_f, cg_args);
  }
};

struct PeriodicSetup{
  typedef ActionsPeriodic ActionsType;
  typedef PeriodicGimplD GimplD;
  typedef CayleyFermion5D<WilsonImplD> ActionTypeD;
  
  inline static std::unique_ptr<ActionsType> setupActions(const MeasArgs &args, LatticeGaugeFieldD &U_d, Grids &gridsD, LatticeGaugeFieldF &U_f, Grids &gridsF){ 
    WilsonImplD::ImplParams action_params = setupActionParamsPeriodic();
    return std::unique_ptr<ActionsPeriodic>(new ActionsPeriodic(args.action, action_params, args.mass, args.mobius_scale, U_d, gridsD, U_f, gridsF));
  }

  inline static void invertRandomSource(FermionFieldPeriodicD &rand_vol_sol, const FermionFieldPeriodicD &rand_vol_src, ActionsPeriodic &actions, const MixedCGargs &cg_args){
    mixedPrecInvertField(rand_vol_sol, rand_vol_src, *actions.action_d, *actions.action_f, cg_args);
  }
};

template<typename Setup>
void run(const MeasArgs &args, int cfg_start, int cfg_step, int cfg_lessthan, bool cps_cfg){
  
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];
  
  Grids gridsD = makeDoublePrecGrids(args.Ls, latt);
  Grids gridsF = makeSinglePrecGrids(args.Ls, latt);

  LatticeGaugeFieldD U_d(gridsD.UGrid);
  LatticeGaugeFieldF U_f(gridsF.UGrid);
  LatticeGaugeFieldD V_d(gridsD.UGrid); //smeared
  
  std::unique_ptr<typename Setup::ActionsType> _actions = Setup::setupActions(args, U_d, gridsD, U_f, gridsF);
  typename Setup::ActionsType &actions = *_actions;

  typedef typename Setup::GimplD GimplD;
  typedef typename Setup::ActionTypeD ActionTypeD;
  typedef typename ActionTypeD::FermionField FermionFieldTypeD;
  
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
      readCPSconfiguration<GimplD>(U_d, sRNG, pRNG, traj, args.cfg_stub) :
      readConfiguration<GimplD>(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);

    precisionChange(U_f,U_d);
    actions.ImportGauge(U_d,U_f);

    //Timeslice plaquette
    {
      auto tslice_plaq = WilsonLoops<GimplD>::timesliceAvgSpatialPlaquette(U_d);
      asciiWriteArray(tslice_plaq, "timeslice_plaq", traj);   
    }    
    //Oriented plaquette
    {
      RealD plaq;
      plaq = avgOrientedPlaquette<GimplD>(0,1,U_d); asciiWriteValue(plaq, "plaq_XY", traj);
      plaq = avgOrientedPlaquette<GimplD>(0,2,U_d); asciiWriteValue(plaq, "plaq_XZ", traj);
      plaq = avgOrientedPlaquette<GimplD>(1,2,U_d); asciiWriteValue(plaq, "plaq_YZ", traj);
      plaq = avgOrientedPlaquette<GimplD>(0,3,U_d); asciiWriteValue(plaq, "plaq_XT", traj);
      plaq = avgOrientedPlaquette<GimplD>(1,3,U_d); asciiWriteValue(plaq, "plaq_YT", traj);
      plaq = avgOrientedPlaquette<GimplD>(2,3,U_d); asciiWriteValue(plaq, "plaq_ZT", traj);
    }

    //Chiral condensate
    {
      FermionFieldTypeD rand_vol_src = randomGaussianVolumeSource<ActionTypeD>(pRNG, gridsD.UGrid);
      FermionFieldTypeD rand_vol_sol(gridsD.UGrid);
      Setup::invertRandomSource(rand_vol_sol, rand_vol_src, actions, args.cg_args);
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
    WilsonFlowAdaptiveMeasGeneral<GimplD>(wflow_io, args.wflow_init_epsilon, args.wflow_maxTau, args.wflow_tolerance, U_d, &V_d);
	
    asciiWriteArray(wflow_io.energy_density_clover, "wflow_clover", traj);
    asciiWriteArray(wflow_io.energy_density_plaq, "wflow_plaq", traj);
    writeTsliceTopQsmr(wflow_io.timeslice_topq, "timeslice_topq5li_smr", traj);
    writeTsliceTopQsmr(wflow_io.timeslice_plaq, "timeslice_plaq_smr", traj);

    //Measure topological charge
    std::vector<std::vector<Real> > tslice_topq5li_contribs = timesliceTopologicalCharge5LiContributions<GimplD>(V_d);
    asciiWriteArray(tslice_topq5li_contribs, "timeslice_topq5li_contribs", traj);

    std::vector<Real> tslice_topq5li = timesliceTopologicalCharge5Li(tslice_topq5li_contribs);
    asciiWriteArray(tslice_topq5li, "timeslice_topq5li", traj);

    std::vector<Real> topq5li_contribs = topologicalCharge5LiContributions(tslice_topq5li_contribs);
    asciiWriteArray(topq5li_contribs, "topq5li_contribs", traj);
   
    RealD topq5li = topologicalCharge5Li(topq5li_contribs);
    asciiWriteValue(topq5li, "topq5li", traj);
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

  int ngp = 0;
  for(int i=0;i<4;i++) ngp += args.GparityDirs[i];

  if(ngp == 0){
    run<PeriodicSetup>(args, cfg_start, cfg_step, cfg_lessthan, cps_cfg);
  }else{
    run<GparitySetup>(args, cfg_start, cfg_step, cfg_lessthan, cps_cfg);
  }
  
  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
