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
				  LanczosParameters, lanc_args_s,
				  int, evec_prec
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
    evec_prec = 2; //2=double 1=single
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


struct EvecContainerOpts{
  bool load_evecs;
  std::string load_evecs_stub;
  std::string load_evals_stub;

  bool save_evecs;
  std::string save_evecs_stub;
  std::string save_evals_stub;

  EvecContainerOpts(): load_evecs(false), save_evecs(false){}
};

struct EvecContainer{
  std::vector<RealD> * eval_ptr;
  std::vector<FermionFieldD> * evec_ptr;
  std::vector<FermionFieldF> * evec_ptr_f;
  int precision;

  EvecContainer(): eval_ptr(nullptr), evec_ptr(nullptr), evec_ptr_f(nullptr), precision(0){}

  void clear(){
    if(eval_ptr){ delete eval_ptr; eval_ptr = nullptr; }
    if(evec_ptr){ delete evec_ptr; evec_ptr = nullptr; }
    if(evec_ptr_f){ delete evec_ptr_f; evec_ptr_f = nullptr;} 
  }
  
  ~EvecContainer(){
    clear();
  }

  //Double prec
  void generate(const LanczosParameters &lanc_arg, int traj, 
		Grids &grids, CayleyFermion5D<GparityWilsonImplD> &action, 
		const LatticeGaugeFieldD &U_d, GridParallelRNG &pRNG, 
		const EvecContainerOpts &opts= EvecContainerOpts() ){
    clear();
    eval_ptr = new std::vector<RealD>();
    evec_ptr = new std::vector<FermionFieldD>();
    precision = 2;

    if(opts.load_evecs){
      readEigenvalues(*eval_ptr, *evec_ptr, grids.FrbGrid, opts.load_evals_stub, opts.load_evecs_stub, traj);
    }else{
      computeEigenvalues<CayleyFermion5D<GparityWilsonImplD>, FermionFieldD>(*eval_ptr, *evec_ptr, lanc_arg, grids.FGrid, grids.FrbGrid, U_d, action, pRNG);
    }
    if(opts.save_evecs) saveEigenvalues(*eval_ptr, *evec_ptr, opts.save_evals_stub, opts.save_evecs_stub, traj);
  }

  //Single
  void generate(const LanczosParameters &lanc_arg, int traj, 
		Grids &grids, CayleyFermion5D<GparityWilsonImplF> &action, 
		const LatticeGaugeFieldF &U_f, GridParallelRNG &pRNG, 
		const EvecContainerOpts &opts= EvecContainerOpts() ){
    clear();
    eval_ptr = new std::vector<RealD>();
    evec_ptr_f = new std::vector<FermionFieldF>();
    precision = 1;

    if(opts.load_evecs){
      readEigenvalues(*eval_ptr, *evec_ptr_f, grids.FrbGrid, opts.load_evals_stub, opts.load_evecs_stub, traj);
    }else{
      computeEigenvalues<CayleyFermion5D<GparityWilsonImplF>, FermionFieldF>(*eval_ptr, *evec_ptr_f, lanc_arg, grids.FGrid, grids.FrbGrid, U_f, action, pRNG);
    }
    if(opts.save_evecs) saveEigenvalues(*eval_ptr, *evec_ptr_f, opts.save_evals_stub, opts.save_evecs_stub, traj);
  }

  void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				  const LatticeSCFmatrixD &msrc, CayleyFermion5D<GparityWilsonImplD> &action_d, CayleyFermion5D<GparityWilsonImplF> &action_f, 
				  double tol, double inner_tol){
    if(precision == 2){
      std::cout << GridLogMessage << "Inverting using double precision evecs" << std::endl;
      assert(eval_ptr != nullptr && evec_ptr != nullptr);
      GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, tol, inner_tol, eval_ptr, evec_ptr);
    }
    else if(precision == 1){
      std::cout << GridLogMessage << "Inverting using single precision evecs" << std::endl;
      assert(eval_ptr != nullptr && evec_ptr_f != nullptr);
      GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, tol, inner_tol, eval_ptr, evec_ptr_f);
    }
    else if(precision == 0){ //no evecs
      std::cout << GridLogMessage << "Inverting without evecs" << std::endl;
      GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, tol, inner_tol);
    }
    else assert(0);
  }

  void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					   const LatticeSCFmatrixD &msrc,
					   CayleyFermion5D<GparityWilsonImplD> &action_d,
					   CayleyFermion5D<GparityWilsonImplD> &subgrid_action_d, CayleyFermion5D<GparityWilsonImplF> &subgrid_action_f,
					   double tol, double inner_tol){
    if(precision == 2){
      std::cout << GridLogMessage << "Inverting using double precision evecs" << std::endl;
      assert(eval_ptr != nullptr && evec_ptr != nullptr);
      GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, eval_ptr, evec_ptr);
    }
    else if(precision == 1){
      std::cout << GridLogMessage << "Inverting using single precision evecs" << std::endl;
      assert(eval_ptr != nullptr && evec_ptr_f != nullptr);
      GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, eval_ptr, evec_ptr_f);
    }
    else if(precision == 0){ //no evecs
      std::cout << GridLogMessage << "Inverting without evecs" << std::endl;
      GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol);
    }
    else assert(0);
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

  EvecContainerOpts evec_opts_l, evec_opts_s;
 
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
      evec_opts_l.save_evecs = true;
      evec_opts_l.save_evecs_stub = argv[i+1];
      evec_opts_l.save_evals_stub = argv[i+2];
    }else if(sargv == "--load_evecs"){
      evec_opts_l.load_evecs = true;
      evec_opts_l.load_evecs_stub = argv[i+1];
      evec_opts_l.load_evals_stub = argv[i+2];
    }else if(sargv == "--save_evecs_strange"){
      evec_opts_s.save_evecs = true;
      evec_opts_s.save_evecs_stub = argv[i+1];
      evec_opts_s.save_evals_stub = argv[i+2];
    }else if(sargv == "--load_evecs_strange"){
      evec_opts_s.load_evecs = true;
      evec_opts_s.load_evecs_stub = argv[i+1];
      evec_opts_s.load_evals_stub = argv[i+2];
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

  assert(args.evec_prec == 2 || args.evec_prec == 1);

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
  std::vector<int> src_t(args.nsrc);
  if(args.nsrc != 0){ //allow 0 as we might want to just generate and save the eigenvectors / gfix config  
    assert(Lt % args.nsrc == 0);
    int tsep = Lt / args.nsrc;
  
    for(int i=0;i<args.nsrc;i++){    
      int t = i * tsep;
      src_t[i] = t;
    }
  }

  //Start calculation
  EvecContainer eval, eval_s;

  printMem("Pre trajectory loop");
  
  //Start traj loop
  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    std::cout << GridLogMessage << "Starting traj " << traj << std::endl;
    std::vector<int> seeds4({traj, traj+2, traj+3, traj+4});

    GridParallelRNG pRNG(GridsD.UGrid); //4D!
    pRNG.SeedFixedIntegers(seeds4);

    GridParallelRNG pRNG_f(GridsF.UGrid);
    pRNG_f.SeedFixedIntegers(seeds4);
    
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
    if(use_evecs){
      std::cout << GridLogMessage << "Obtaining light eigenvectors" << std::endl;
      if(args.evec_prec == 2)
	eval.generate(args.lanc_args, traj, GridsD, *actions.light.action_d, U_d, pRNG, evec_opts_l);
      else if(args.evec_prec == 1)
	eval.generate(args.lanc_args, traj, GridsF, *actions.light.action_f, U_f, pRNG_f, evec_opts_l);
      else assert(0);
    }else{
      std::cout << GridLogMessage << "Not using light eigenvectors" << std::endl;
    }

    printMem("Post light evec generation");

    //////////////////////////////////// Strange evecs ///////////////////////////////////////////
    if(use_evecs_strange){
      std::cout << GridLogMessage << "Obtaining strange eigenvectors" << std::endl;
      if(args.evec_prec == 2)
	eval.generate(args.lanc_args, traj, GridsD, *actions.strange.action_d, U_d, pRNG, evec_opts_s);
      else if(args.evec_prec == 1)
	eval.generate(args.lanc_args, traj, GridsF, *actions.strange.action_f, U_f, pRNG_f, evec_opts_s);
      else assert(0);
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
      if(use_split_grid)
	eval.splitGridMixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d,  *actions_sub.light.action_d, *actions_sub.light.action_f, args.cg_stop, args.cg_stop_inner);
      else
	eval.mixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d, *actions.light.action_f, args.cg_stop, args.cg_stop_inner);
      
      const LatticeSCFmatrixD &Rp2 = Rp1;
      const LatticeSCFmatrixD &Rp2_mid = Rp1_mid;

      //Do strange quark also
      std::cout << GridLogMessage << "Starting strange quark inverse" << std::endl;
      LatticeSCFmatrixD Rp1_s(GridsD.UGrid), Rp1_s_mid(GridsD.UGrid);
      if(use_split_grid)
	eval_s.splitGridMixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d,  *actions_sub.strange.action_d, *actions_sub.strange.action_f, args.cg_stop, args.cg_stop_inner);
      else
	eval_s.mixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d, *actions.strange.action_f, args.cg_stop, args.cg_stop_inner);      

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

    if(args.nsrc != 0){
      asciiWriteArray(Ct_pion, "pion_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_j5q, "j5q_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_ps_singlet, "ps_singlet_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_kaon, "kaon_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_j5q_kaon, "j5q_kaon_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
    }
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
