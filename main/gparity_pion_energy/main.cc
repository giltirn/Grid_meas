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
				  MixedCGargs, cg_args,
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
    cg_args.algorithm = MixedCGalgorithm::RestartedCG;
    cg_args.tolerance = 1e-8;
    cg_args.restartedcg_inner_tol = 1e-6;
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

template<bool isGparity>
struct _action_call_wrapper{};

template<>
struct _action_call_wrapper<true>{
  template<typename ActionD, typename ActionF, typename EvecFieldType>
  static void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
					 const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
					 const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
    std::cout << GridLogMessage << "Inverting with G-parity action" << std::endl;
    std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
    std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
    GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
  }
  template<typename ActionD, typename ActionF, typename EvecFieldType>
  static void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						  const LatticeSCFmatrixD &msrc,
						  ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
						  const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
    std::cout << GridLogMessage << "Inverting with G-parity action" << std::endl;
    std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
    std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
    GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval_ptr, evec_ptr);
  }
};
template<>
struct _action_call_wrapper<false>{
  template<typename ActionD, typename ActionF, typename EvecFieldType>
  static void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
					 const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
					 const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
    std::cout << GridLogMessage <<"Inverting with X-conjugate action" << std::endl;
    std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
    std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
    GridMeas::mixedPrecInvertWithMidPropXconj(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
  }
  template<typename ActionD, typename ActionF, typename EvecFieldType>
  static void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						  const LatticeSCFmatrixD &msrc,
						  ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
						  const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
    std::cout << GridLogMessage << "Inverting with X-conjugate action" << std::endl;
    assert(0);
  }
};


template<typename _LanczosAction>
struct EvecContainer{
  typedef _LanczosAction LanczosAction;
  typedef typename LanczosAction::FermionField EvecFieldType;
  
  std::vector<RealD> eval;
  std::vector<EvecFieldType> evec;

  EvecContainer(){}

  void clear(){
    eval.clear();
    std::vector<EvecFieldType>().swap(evec);
  }

  void generate(const LanczosParameters &lanc_arg, int traj, 
		LanczosAction &action, 
		const typename LanczosAction::GaugeField &U, GridParallelRNG &pRNG, 
		const EvecContainerOpts &opts= EvecContainerOpts() ){
    GridCartesian* FGrid = (GridCartesian*)action.FermionGrid();
    GridRedBlackCartesian* FrbGrid = (GridRedBlackCartesian*)action.FermionRedBlackGrid();

    clear();
    if(opts.load_evecs){
      readEigenvalues(eval, evec, FrbGrid, opts.load_evals_stub, opts.load_evecs_stub, traj);
    }else{
      computeEigenvalues<LanczosAction, EvecFieldType>(eval, evec, lanc_arg, FGrid, FrbGrid, U, action, pRNG);
    }
    if(opts.save_evecs) saveEigenvalues(eval, evec, opts.save_evals_stub, opts.save_evecs_stub, traj);
  }
  
  template<typename ActionD, typename ActionF>
  void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				  const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
				  const MixedCGargs &args) const{
    _action_call_wrapper<ActionD::isGparity>::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, eval, evec);
  }
  template<typename ActionD, typename ActionF>
  void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					   const LatticeSCFmatrixD &msrc,
					   ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
					   const MixedCGargs &args) const{
    _action_call_wrapper<ActionD::isGparity>::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval, evec);
  }
};

struct Opts{
  std::string arg_file;
  int cfg_start;
  int cfg_lessthan;
  int cfg_step;

  EvecContainerOpts evec_opts_l;
  EvecContainerOpts evec_opts_s;
 
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

  Opts(){
    load_gfix_cfg = false;
    save_gfix_cfg = false;
    cps_cfg = false;
    unit_gauge = false;
    use_evecs = true;
    use_evecs_strange = true;

    use_split_grid = false;
    use_gauge_fixing = true;
  }
};

inline void addResult(std::vector<RealD> &into, const std::vector<RealD> &from, const int nsrc){
  for(int t=0;t<into.size();t++) //average over sources
    into[t] += from[t] / RealD(nsrc);
}

template<typename LanczosAction>
void run(const MeasArgs &args, const Opts &opts){
  printMem("Program body start");
  
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];

  std::cout << GridLogMessage << "Creating 4d and 5d Grids with Ls=" << args.Ls << std::endl;
  Grids GridsD = makeDoublePrecGrids(args.Ls, latt);
  Grids GridsF = makeSinglePrecGrids(args.Ls, latt);

  printMem("Post Grid creation");
  
  GparityWilsonImplD::ImplParams Params = setupActionParams(args.GparityDirs);

  LatticeGaugeFieldD U_d(GridsD.UGrid);
  LatticeGaugeFieldF U_f(GridsF.UGrid);
  typename LanczosAction::GaugeField & U_lanczos = getInstance<typename LanczosAction::GaugeField>(U_d,U_f);

  printMem("Post gauge field creation");
  
  //Get the actions
  ActionsLightStrange actions(args.action, args.ml, args.ms, Params, args.mobius_scale, U_d, GridsD, U_f, GridsF);

  printMem("Post action creation");
  
  //Prepare split grid if in use
  Grids *SubGridsD_p = nullptr, *SubGridsF_p = nullptr;
  ActionsLightStrange actions_sub;
  LatticeGaugeFieldD *U_sub_d_p = nullptr;
  LatticeGaugeFieldF *U_sub_f_p = nullptr;
  if(opts.use_split_grid){
    printMem("Pre split Grid creation");
    Coordinate proc_sub(opts.split_grid_proc_layout);
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

  printMem("Pre trajectory loop");
  
  //Start traj loop
  for(int traj = opts.cfg_start; traj < opts.cfg_lessthan; traj += opts.cfg_step){
    std::cout << GridLogMessage << "Starting traj " << traj << std::endl;
    std::vector<int> seeds4({traj, traj+2, traj+3, traj+4});

    GridParallelRNG pRNG(GridsD.UGrid); //4D!
    pRNG.SeedFixedIntegers(seeds4);

    GridParallelRNG pRNG_f(GridsF.UGrid);
    pRNG_f.SeedFixedIntegers(seeds4);
    
    GridParallelRNG &pRNG_lanczos = args.evec_prec == 2 ? pRNG : pRNG_f;

    GridSerialRNG sRNG;  
    sRNG.SeedFixedIntegers(seeds4); 
    
    ///////////////////////////////// Read config /////////////////////////////////////////
    if(opts.unit_gauge){
      std::cout << GridLogMessage << "Setting gauge field to unit gauge" << std::endl;
      SU<Nc>::ColdConfiguration(U_d);
    }else{
      std::cout << GridLogMessage << "Reading configuration" << std::endl;
      opts.cps_cfg ? 
	readCPSconfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub) :
	readConfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);
    }
    printMem("Post configuration read");
    
    ///////////////////////////////// Gauge fix ///////////////////////////////////////////
    if(opts.use_gauge_fixing){
      if(opts.load_gfix_cfg)
	readGaugeFixedConfiguration(U_d, opts.load_gfix_cfg_stub, traj);    
      else
	CoulombGaugeFix(U_d, args.gfix_alpha, args.gfix_stop_cnd, args.gfix_fourier_accelerate);
      
      if(opts.save_gfix_cfg)
	writeGaugeFixedConfiguration(U_d, opts.save_gfix_cfg_stub, traj);
    }
    printMem("Post gauge fix");
    
    precisionChange(U_f, U_d);

    actions.ImportGauge(U_d,U_f);

    if(opts.use_split_grid){
      Grid_split(U_d,*U_sub_d_p);
      precisionChange(*U_sub_f_p,*U_sub_d_p);
      actions_sub.ImportGauge(*U_sub_d_p,*U_sub_f_p);      
    }
    printMem("Post gauge field import");
    
    //Start calculation
    //////////////////////////////////// Light evecs ///////////////////////////////////////////
    typedef EvecContainer<LanczosAction> EvecContainerType;
    EvecContainerType eval;
    if(opts.use_evecs){
      std::cout << GridLogMessage << "Obtaining light eigenvectors" << std::endl;
      eval.generate(args.lanc_args, traj, *actions.light.getAction<LanczosAction>(), U_lanczos, pRNG_lanczos, opts.evec_opts_l);
    }else{
      std::cout << GridLogMessage << "Not using light eigenvectors" << std::endl;
    }

    printMem("Post light evec generation");

    //////////////////////////////////// Strange evecs ///////////////////////////////////////////
    EvecContainerType eval_s;
    if(opts.use_evecs_strange){
      std::cout << GridLogMessage << "Obtaining strange eigenvectors" << std::endl;
      eval_s.generate(args.lanc_args_s, traj, *actions.strange.getAction<LanczosAction>(), U_lanczos, pRNG_lanczos, opts.evec_opts_s);
    }else{
      std::cout << GridLogMessage << "Not using strange eigenvectors" << std::endl;
    }

    printMem("Post strange evec generation");
    
    /////////////////////////////////// Trajectory loop //////////////////////////////////////////
    std::vector<RealD> Ct_pion_momwall(Lt, 0);
    std::vector<RealD> Ct_j5q_momwall(Lt, 0);
    std::vector<RealD> Ct_kaon_momwall(Lt, 0);
    std::vector<RealD> Ct_ps_singlet_momwall(Lt, 0);
    std::vector<RealD> Ct_j5q_kaon_momwall(Lt, 0);

    std::vector<RealD> Ct_pion_coswall(Lt, 0);
    std::vector<RealD> Ct_j5q_coswall(Lt, 0);
    std::vector<RealD> Ct_kaon_coswall(Lt, 0);
    std::vector<RealD> Ct_ps_singlet_coswall(Lt, 0);
    std::vector<RealD> Ct_j5q_kaon_coswall(Lt, 0);

    for(int s=0;s<args.nsrc;s++){ 
      int t0 = src_t[s];
      std::cout << GridLogMessage << "Starting calculation with source timeslice t0=" << t0 << std::endl;

      {
	//Coulomb gauge fixed wall momentum sources (must be inverted with regular G-parity Dirac op)
	LatticeSCFmatrixD src_p1 = momentumWallSource(p1_phys, t0, GridsD.UGrid);
	LatticeSCFmatrixD Rp1(GridsD.UGrid), Rp1_mid(GridsD.UGrid);
	std::cout << GridLogMessage << "Starting light quark inverse" << std::endl;
	if(opts.use_split_grid)
	  eval.splitGridMixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d, *actions_sub.light.action_d, 
						   *actions_sub.light.action_f, args.cg_args);
	else
	  eval.mixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.action_d, *actions.light.action_f, args.cg_args);
      
	const LatticeSCFmatrixD &Rp2 = Rp1;
	const LatticeSCFmatrixD &Rp2_mid = Rp1_mid;

	//Do strange quark also
	std::cout << GridLogMessage << "Starting strange quark inverse" << std::endl;
	LatticeSCFmatrixD Rp1_s(GridsD.UGrid), Rp1_s_mid(GridsD.UGrid);
	if(opts.use_split_grid)
	  eval_s.splitGridMixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d, *actions_sub.strange.action_d, 
						     *actions_sub.strange.action_f, args.cg_args);
	else
	  eval_s.mixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.action_d, *actions.strange.action_f, args.cg_args);      

	addResult(Ct_pion_momwall, momWallSourcePionCorrelator(p1, p2, t0, Rp1, Rp2), args.nsrc);
	addResult(Ct_j5q_momwall, momWallSourcePionCorrelator(p1, p2, t0, Rp1_mid, Rp2_mid), args.nsrc);
	addResult(Ct_ps_singlet_momwall, momWallSourcePSsingletCorrelator(p1, t0, Rp1), args.nsrc);
	addResult(Ct_kaon_momwall, momWallSourceKaonCorrelator(p1, t0, Rp1, Rp1_s), args.nsrc);
	addResult(Ct_j5q_kaon_momwall, momWallSourceKaonCorrelator(p1, t0, Rp1_mid, Rp1_s_mid), args.nsrc);
      }


      {
	//Coulomb gauge fixed cosine wall momentum sources, use X-conjugate Dirac operator
	LatticeSCFmatrixD src_p1 = cosineWallSource(p1_phys, t0, GridsD.UGrid);
	LatticeSCFmatrixD Rp1(GridsD.UGrid), Rp1_mid(GridsD.UGrid);
	std::cout << GridLogMessage << "Starting light quark inverse" << std::endl;
	if(opts.use_split_grid)
	  eval.splitGridMixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.xconj_action_d, *actions_sub.light.xconj_action_d, 
						   *actions_sub.light.xconj_action_f, args.cg_args);
	else
	  eval.mixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *actions.light.xconj_action_d, *actions.light.xconj_action_f, args.cg_args);
      
	const LatticeSCFmatrixD &Rp2 = Rp1;
	const LatticeSCFmatrixD &Rp2_mid = Rp1_mid;

	//Do strange quark also
	std::cout << GridLogMessage << "Starting strange quark inverse" << std::endl;
	LatticeSCFmatrixD Rp1_s(GridsD.UGrid), Rp1_s_mid(GridsD.UGrid);
	if(opts.use_split_grid)
	  eval_s.splitGridMixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.xconj_action_d, *actions_sub.strange.xconj_action_d, 
						     *actions_sub.strange.xconj_action_f, args.cg_args);
	else
	  eval_s.mixedPrecInvertWithMidProp(Rp1_s, Rp1_s_mid, src_p1, *actions.strange.xconj_action_d, *actions.strange.xconj_action_f, args.cg_args);      

	addResult(Ct_pion_coswall, momWallSourcePionCorrelator(p1, p2, t0, Rp1, Rp2), args.nsrc);
	addResult(Ct_j5q_coswall, momWallSourcePionCorrelator(p1, p2, t0, Rp1_mid, Rp2_mid), args.nsrc);
	addResult(Ct_ps_singlet_coswall, momWallSourcePSsingletCorrelator(p1, t0, Rp1), args.nsrc);
	addResult(Ct_kaon_coswall, momWallSourceKaonCorrelator(p1, t0, Rp1, Rp1_s), args.nsrc);
	addResult(Ct_j5q_kaon_coswall, momWallSourceKaonCorrelator(p1, t0, Rp1_mid, Rp1_s_mid), args.nsrc);
      }

    }

    if(args.nsrc != 0){
      asciiWriteArray(Ct_pion_momwall, "pion_momwall_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_j5q_momwall, "j5q_momwall_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_ps_singlet_momwall, "ps_singlet_momwall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_kaon_momwall, "kaon_momwall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_j5q_kaon_momwall, "j5q_kaon_momwall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);

      asciiWriteArray(Ct_pion_coswall, "pion_coswall_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_j5q_coswall, "j5q_coswall_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
      asciiWriteArray(Ct_ps_singlet_coswall, "ps_singlet_coswall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_kaon_coswall, "kaon_coswall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
      asciiWriteArray(Ct_j5q_kaon_coswall, "j5q_kaon_coswall_mom" + momstr(p1) + "_mom" + momstr(mp1), traj);
    }
  }
}

  
int main(int argc, char** argv){
  Grid_init(&argc, &argv);
  printMem("Program start");
  
  assert(argc >= 5);
  std::string arg_file = argv[1];
  
  Opts opts;
  opts.cfg_start = std::stoi(argv[2]);
  opts.cfg_lessthan = std::stoi(argv[3]);
  opts.cfg_step = std::stoi(argv[4]);
  
  for(int i=1;i<argc;i++){
    std::string sargv(argv[i]);
    if(sargv == "--cps_cfg"){
      opts.cps_cfg = true;
    }else if(sargv == "--save_evecs"){
      opts.evec_opts_l.save_evecs = true;
      opts.evec_opts_l.save_evecs_stub = argv[i+1];
      opts.evec_opts_l.save_evals_stub = argv[i+2];
    }else if(sargv == "--load_evecs"){
      opts.evec_opts_l.load_evecs = true;
      opts.evec_opts_l.load_evecs_stub = argv[i+1];
      opts.evec_opts_l.load_evals_stub = argv[i+2];
    }else if(sargv == "--save_evecs_strange"){
      opts.evec_opts_s.save_evecs = true;
      opts.evec_opts_s.save_evecs_stub = argv[i+1];
      opts.evec_opts_s.save_evals_stub = argv[i+2];
    }else if(sargv == "--load_evecs_strange"){
      opts.evec_opts_s.load_evecs = true;
      opts.evec_opts_s.load_evecs_stub = argv[i+1];
      opts.evec_opts_s.load_evals_stub = argv[i+2];
    }else if(sargv == "--load_gfix_cfg"){
      opts.load_gfix_cfg=true;
      opts.load_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--save_gfix_cfg"){
      opts.save_gfix_cfg=true;
      opts.save_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--disable_evecs"){
      opts.use_evecs = false;
    }else if(sargv == "--disable_evecs_strange"){
      opts.use_evecs_strange = false;
    }else if(sargv == "--unit_gauge"){
      opts.unit_gauge = true;
    }else if(sargv == "--split_grid"){
      opts.use_split_grid = true;
      GridCmdOptionIntVector(argv[i+1], opts.split_grid_proc_layout);
      assert(opts.split_grid_proc_layout.size() == Nd);
      std::cout << GridLogMessage << "Using split grid with geometry (";
      for(int d=0;d<Nd;d++) std::cout << opts.split_grid_proc_layout[d] << " ";
      std::cout << ")" << std::endl;
    }else if(sargv == "--disable_gauge_fixing"){
      opts.use_gauge_fixing = false;
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

  if(args.evec_prec == 2){
    run<CayleyFermion5D<XconjugateWilsonImplD> >(args,opts);   
    //run<CayleyFermion5D<GparityWilsonImplD> >(args,opts);   
  }else if(args.evec_prec == 1){
    run<CayleyFermion5D<XconjugateWilsonImplF> >(args,opts);
    //run<CayleyFermion5D<GparityWilsonImplF> >(args,opts);   
  }else{
    assert(0);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
