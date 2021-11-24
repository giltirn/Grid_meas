#include<Grid/Grid.h>
#include <common.h>

using namespace GridMeas;
using namespace Grid;

struct MeasArgs: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MeasArgs,
				  ActionType, action,
				  int, Ls,
				  double, mobius_scale,
				  std::vector<Integer>, GparityDirs,
				  std::string, cfg_stub,
				  std::string, rng_stub,
				  double, gfix_alpha,
				  double, gfix_stop_cnd,
				  bool, gfix_fourier_accelerate,
				  double, ml,
				  double, ms,
				  double, cg_stop,
				  double, cg_stop_inner,
				  int, nsrc,
				  LanczosParameters, lanc_args,
				  LanczosParameters, lanc_args_s
				  );
  MeasArgs() { 
    action = ActionType::DWF;
    Ls = 12;
    GparityDirs = {1,0,0};
    cfg_stub = "ckpoint_lat";
    rng_stub = "ckpoint_rng";
    gfix_alpha = 0.05;
    gfix_stop_cnd = 1e-9;
    gfix_fourier_accelerate = false;    
    ml = 0.01;
    ms = 0.032;
    mobius_scale = 2.0;
    cg_stop = 1e-8;
    cg_stop_inner = 1e-6;
    nsrc = 32;
  }
};

struct Actions{
  CayleyFermion5D<GparityWilsonImplD>* action_d;
  CayleyFermion5D<GparityWilsonImplF>* action_f;
  
  Actions(): action_d(nullptr), action_f(nullptr){}

  Actions(ActionType action,
	  const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
	  LatticeGaugeFieldD &U_d, Grids &GridsD,
	  LatticeGaugeFieldF &U_f, Grids &GridsF){
    printMem("Pre actionD creation");
    action_d = createActionD(action, Params, mass, mobius_scale, U_d, GridsD);
    printMem("Pre actionF creation");
    action_f = createActionF(action, Params, mass, mobius_scale, U_f, GridsF);
    printMem("Post actionD/F creation");
  }    
  
  ~Actions(){
    if(action_d) delete action_d;
    if(action_f) delete action_f;
  }

  void ImportGauge(LatticeGaugeFieldD &U_d,LatticeGaugeFieldF &U_f){
    if(action_d) action_d->ImportGauge(U_d);
    if(action_f) action_f->ImportGauge(U_f);
  }
		   
    
};

struct ActionsLightStrange{
  Actions light;
  Actions strange;

  ActionsLightStrange(){}
  
  ActionsLightStrange(ActionType action,
		      double light_mass, double strange_mass,
		      const GparityWilsonImplD::ImplParams &Params, double mobius_scale,
		      LatticeGaugeFieldD &U_d, Grids &GridsD,
		      LatticeGaugeFieldF &U_f, Grids &GridsF):
    light(action,Params,light_mass,mobius_scale,U_d,GridsD,U_f,GridsF),
    strange(action,Params,strange_mass,mobius_scale,U_d,GridsD,U_f,GridsF)
  {}

  void ImportGauge(LatticeGaugeFieldD &U_d,LatticeGaugeFieldF &U_f){
    light.ImportGauge(U_d,U_f);
    strange.ImportGauge(U_d,U_f);
  }

};


  
int main(int argc, char** argv){
  Grid_init(&argc, &argv);
  printMem("Program start");
  
  assert(argc >= 5);
  std::string arg_file = argv[1];
  int cfg_start = std::stoi(argv[2]);
  int cfg_lessthan = std::stoi(argv[3]);
  int cfg_step = std::stoi(argv[4]);
 
  bool load_evecs = false;
  std::string load_evecs_stub, load_evals_stub;

  bool save_evecs = false;
  std::string save_evecs_stub, save_evals_stub;

  bool load_evecs_strange = false;
  std::string load_evecs_strange_stub, load_evals_strange_stub;

  bool save_evecs_strange = false;
  std::string save_evecs_strange_stub, save_evals_strange_stub;

  bool load_gfix_cfg = false;
  std::string load_gfix_cfg_stub;

  bool save_gfix_cfg = false;
  std::string save_gfix_cfg_stub;

  bool cps_cfg = false;

  bool unit_gauge = false;
  bool use_evecs = true;
  bool use_evecs_strange = true;

  bool use_split_grid = false;
  std::vector<int> split_grid_proc_layout;

  bool use_gauge_fixing = true;
  
  for(int i=1;i<argc;i++){
    std::string sargv(argv[i]);
    if(sargv == "--cps_cfg"){
      cps_cfg = true;
    }else if(sargv == "--save_evecs"){
      save_evecs = true;
      save_evecs_stub = argv[i+1];
      save_evals_stub = argv[i+2];
    }else if(sargv == "--load_evecs"){
      load_evecs = true;
      load_evecs_stub = argv[i+1];
      load_evals_stub = argv[i+2];
    }else if(sargv == "--save_evecs_strange"){
      save_evecs_strange = true;
      save_evecs_strange_stub = argv[i+1];
      save_evals_strange_stub = argv[i+2];
    }else if(sargv == "--load_evecs_strange"){
      load_evecs_strange = true;
      load_evecs_strange_stub = argv[i+1];
      load_evals_strange_stub = argv[i+2];
    }else if(sargv == "--load_gfix_cfg"){
      load_gfix_cfg=true;
      load_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--save_gfix_cfg"){
      save_gfix_cfg=true;
      save_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--disable_evecs"){
      use_evecs = false;
    }else if(sargv == "--disable_evecs_strange"){
      use_evecs_strange = false;
    }else if(sargv == "--unit_gauge"){
      unit_gauge = true;
    }else if(sargv == "--split_grid"){
      use_split_grid = true;
      GridCmdOptionIntVector(argv[i+1], split_grid_proc_layout);
      assert(split_grid_proc_layout.size() == Nd);
      std::cout << GridLogMessage << "Using split grid with geometry (";
      for(int d=0;d<Nd;d++) std::cout << split_grid_proc_layout[d] << " ";
      std::cout << ")" << std::endl;
    }else if(sargv == "--disable_gauge_fixing"){
      use_gauge_fixing = false;
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
    {
      Grid::XmlWriter wr(arg_file + ".templ");
      write(wr, "Args", args);
    }
      
    std::cout << GridLogMessage << " Done" << std::endl;
    Grid_finalize();
    return 0;
  }

  printMem("Program body start");
  
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];

  std::cout << GridLogMessage << "Creating 4d and 5d Grids with Ls=" << args.Ls << std::endl;
  Grids GridsD = makeDoublePrecGrids(args.Ls, latt);
  Grids GridsF = makeSinglePrecGrids(args.Ls, latt);

  printMem("Post Grid creation");
  
  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = args.GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = args.GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction
  std::cout << GridLogMessage << "Fermion BCs: " << Params.twists << std::endl;

  LatticeGaugeFieldD U_d(GridsD.UGrid);
  LatticeGaugeFieldF U_f(GridsF.UGrid);

  printMem("Post gauge field creation");
  
  ActionsLightStrange actions(args.action, args.ml, args.ms, Params, args.mobius_scale, U_d, GridsD, U_f, GridsF);

  printMem("Post action creation");
  
  //Prepare split grid if in use
  Grids *SubGridsD_p = nullptr, *SubGridsF_p = nullptr;
  ActionsLightStrange actions_sub;
  LatticeGaugeFieldD *U_sub_d_p = nullptr;
  LatticeGaugeFieldF *U_sub_f_p = nullptr;
  if(use_split_grid){
    printMem("Pre split Grid creation");
    Coordinate proc_sub(split_grid_proc_layout);
    SubGridsD_p = new Grids(makeSplitGrids(GridsD, proc_sub));
    SubGridsF_p = new Grids(makeSplitGrids(GridsF, proc_sub));
    printMem("Post split Grid creation");
    U_sub_d_p = new LatticeGaugeFieldD(SubGridsD_p->UGrid);
    U_sub_f_p = new LatticeGaugeFieldF(SubGridsF_p->UGrid);
    printMem("Post split Grid gauge field creation");
    actions_sub = ActionsLightStrange(args.action, args.ml, args.ms, Params, args.mobius_scale, *U_sub_d_p, *SubGridsD_p, *U_sub_f_p, *SubGridsF_p);
    printMem("Post split-grid Grid/action creation");
  }
  
  //Not going to worry about a rotationally invariant operator,
  //just use p_j = pi/2L for both quark and antiquark in direction j
  std::vector<int> p1(3,0);
  for(int i=0;i<3;i++)
    p1[i] = args.GparityDirs[i] ? 1 : 0;
  std::vector<double> p1_phys = getPhysicalMomentum(p1);
  
  std::vector<int> p2 = p1;
  std::vector<double> p2_phys = p1_phys;

  std::vector<int> mp1 = {-p1[0],-p1[1],-p1[2]};

  //Phase factors
  LatticeComplexD p1_src_phase_field = phaseField(p1_phys, GridsD.UGrid); //exp(-i \vec p \cdot \vec x)
  p1_src_phase_field = conjugate(p1_src_phase_field); //exp(+i \vec p \cdot \vec x)

  const LatticeComplexD &p2_src_phase_field = p1_src_phase_field;

  printMem("Post phase factor generation");
  
  //Wall source timeslices
  assert(Lt % args.nsrc == 0);
  int tsep = Lt / args.nsrc;
  std::vector<int> src_t(args.nsrc);
  for(int i=0;i<args.nsrc;i++){    
    int t = i * tsep;
    src_t[i] = t;
  }

  //Start calculation
  std::vector<RealD> eval, eval_s;
  std::vector<FermionFieldD> evec, evec_s;

  printMem("Pre trajectory loop");
  
  //Start traj loop
  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    std::cout << GridLogMessage << "Starting traj " << traj << std::endl;
    std::vector<int> seeds4({traj, traj+2, traj+3, traj+4});
    GridParallelRNG pRNG(GridsD.UGrid); //4D!
    pRNG.SeedFixedIntegers(seeds4);
    
    GridSerialRNG sRNG;  
    sRNG.SeedFixedIntegers(seeds4); 
    
    ///////////////////////////////// Read config /////////////////////////////////////////
    if(unit_gauge){
      std::cout << GridLogMessage << "Setting gauge field to unit gauge" << std::endl;
      SU<Nc>::ColdConfiguration(U_d);
    }else{
      std::cout << GridLogMessage << "Reading configuration" << std::endl;
      cps_cfg ? 
	readCPSconfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub) :
	readConfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);
    }
    printMem("Post configuration read");
    
    ///////////////////////////////// Gauge fix ///////////////////////////////////////////
    if(use_gauge_fixing){
      if(load_gfix_cfg)
	readGaugeFixedConfiguration(U_d, load_gfix_cfg_stub, traj);    
      else
	CoulombGaugeFix(U_d, args.gfix_alpha, args.gfix_stop_cnd, args.gfix_fourier_accelerate);
      
      if(save_gfix_cfg)
	writeGaugeFixedConfiguration(U_d, save_gfix_cfg_stub, traj);
    }
    printMem("Post gauge fix");
    
    precisionChange(U_f, U_d);

    actions.ImportGauge(U_d,U_f);

    if(use_split_grid){
      Grid_split(U_d,*U_sub_d_p);
      precisionChange(*U_sub_f_p,*U_sub_d_p);
      actions_sub.ImportGauge(*U_sub_d_p,*U_sub_f_p);      
    }
    printMem("Post gauge field import");
    
    //////////////////////////////////// Light evecs ///////////////////////////////////////////
    std::vector<RealD> const* eval_ptr = nullptr;
    std::vector<FermionFieldD> const* evec_ptr = nullptr;
    if(use_evecs){
      std::cout << GridLogMessage << "Obtaining light eigenvectors" << std::endl;
      if(load_evecs){
	readEigenvalues(eval, evec, GridsD.FrbGrid, load_evals_stub, load_evecs_stub, traj);
      }else{
	computeEigenvalues<CayleyFermion5D<GparityWilsonImplD>, FermionFieldD>(eval, evec, args.lanc_args, GridsD.FGrid, GridsD.FrbGrid, U_d, *actions.light.action_d, pRNG);
      }
      if(save_evecs) saveEigenvalues(eval, evec, save_evals_stub, save_evecs_stub, traj);

      eval_ptr = &eval;
      evec_ptr = &evec;
    }else{
      std::cout << GridLogMessage << "Not using light eigenvectors" << std::endl;
    }
    printMem("Post light evec generation");

    //////////////////////////////////// Strange evecs ///////////////////////////////////////////
    std::vector<RealD> const* eval_s_ptr = nullptr;
    std::vector<FermionFieldD> const* evec_s_ptr = nullptr;
    if(use_evecs_strange){
      std::cout << GridLogMessage << "Obtaining strange eigenvectors" << std::endl;
      if(load_evecs_strange){
	readEigenvalues(eval_s, evec_s, GridsD.FrbGrid, load_evals_strange_stub, load_evecs_strange_stub, traj);
      }else{
	computeEigenvalues<CayleyFermion5D<GparityWilsonImplD>, FermionFieldD>(eval_s, evec_s, args.lanc_args_s, GridsD.FGrid, GridsD.FrbGrid, U_d, *actions.strange.action_d, pRNG);
      }
      if(save_evecs_strange) saveEigenvalues(eval_s, evec_s, save_evals_strange_stub, save_evecs_strange_stub, traj);

      eval_s_ptr = &eval_s;
      evec_s_ptr = &evec_s;

    }else{
      std::cout << GridLogMessage << "Not using strange eigenvectors" << std::endl;
    }
    printMem("Post strange evec generation");
    
    /////////////////////////////////// Trajectory loop //////////////////////////////////////////
    std::vector<RealD> Ct_pion(Lt, 0);
    std::vector<RealD> Ct_j5q(Lt, 0);
    std::vector<RealD> Ct_kaon(Lt, 0);
    std::vector<RealD> Ct_ps_singlet(Lt, 0);
    std::vector<RealD> Ct_j5q_kaon(Lt, 0);

    for(int s=0;s<args.nsrc;s++){ 
      int t0 = src_t[s];
      std::cout << GridLogMessage << "Starting calculation with source timeslice t0=" << t0 << std::endl;

      //Coulomb gauge fixed wall momentum sources
      LatticeSCFmatrixD eta = wallSource(t0, GridsD.UGrid);
      LatticeSCFmatrixD src_p1 = p1_src_phase_field * eta;

      LatticeSCFmatrixD Rp1(GridsD.UGrid), Rp1_mid(GridsD.UGrid);
      std::cout << GridLogMessage << "Starting light quark inverse" << std::endl;
      if(use_split_grid){
	splitGridMixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d,  *actions_sub.light.action_d, *actions_sub.light.action_f, args.cg_stop, args.cg_stop_inner, eval_ptr, evec_ptr);
      }else{      
	mixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d, *actions.light.action_f, args.cg_stop, args.cg_stop_inner, eval_ptr, evec_ptr);
      }

      const LatticeSCFmatrixD &Rp2 = Rp1;
      const LatticeSCFmatrixD &Rp2_mid = Rp1_mid;

      //Do strange quark also
      std::cout << GridLogMessage << "Starting strange quark inverse" << std::endl;
      LatticeSCFmatrixD Rp1_s(GridsD.UGrid), Rp1_s_mid(GridsD.UGrid);
      if(use_split_grid){
	splitGridMixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d,  *actions_sub.strange.action_d, *actions_sub.strange.action_f, args.cg_stop, args.cg_stop_inner, eval_ptr, evec_ptr);
      }else{      
	mixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d, *actions.strange.action_f, args.cg_stop, args.cg_stop_inner, eval_s_ptr, evec_s_ptr);
      }

      std::vector<RealD> Cts_pion = momWallSourcePionCorrelator(p1, p2, t0, Rp1, Rp2);
      std::vector<RealD> Cts_j5q = momWallSourcePionCorrelator(p1, p2, t0, Rp1_mid, Rp2_mid);
      std::vector<RealD> Cts_ps_singlet = momWallSourcePSsingletCorrelator(p1, t0, Rp1);
      std::vector<RealD> Cts_kaon = momWallSourceKaonCorrelator(p1, t0, Rp1, Rp1_s);
      std::vector<RealD> Cts_j5q_kaon = momWallSourceKaonCorrelator(p1, t0, Rp1_mid, Rp1_s_mid);

      for(int t=0;t<Lt;t++){ //average over sources
	Ct_pion[t] += Cts_pion[t] / RealD(args.nsrc);
	Ct_j5q[t] += Cts_j5q[t] / RealD(args.nsrc);
	Ct_kaon[t] += Cts_kaon[t] / RealD(args.nsrc);
	Ct_ps_singlet[t] += Cts_ps_singlet[t] / RealD(args.nsrc);
	Ct_j5q_kaon[t] += Cts_j5q_kaon[t] / RealD(args.nsrc);
      }
      printMem("Post measurement");      
    }

    asciiWriteArray(Ct_pion, "pion_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
    asciiWriteArray(Ct_j5q, "j5q_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
    asciiWriteArray(Ct_ps_singlet, "ps_singlet_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
    asciiWriteArray(Ct_kaon, "kaon_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
    asciiWriteArray(Ct_j5q_kaon, "j5q_kaon_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
