#include<Grid/Grid.h>
#include <common.h>

using namespace GridMeas;
using namespace Grid;

struct MeasArgs: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MeasArgs,
				  ActionType, action,
				  int, Ls,
				  double, mobius_scale,
				  std::string, cfg_stub,
				  std::string, rng_stub,
				  double, gfix_alpha,
				  double, gfix_stop_cnd,
				  bool, gfix_fourier_accelerate,
				  double, mass,
				  MixedCGargs, cg_args_sloppy,
				  MixedCGargs, cg_args_exact,
				  int, nsloppy,
				  int, nexact,
				  LanczosParameters, lanc_args,
				  int, evec_prec,
				  std::vector<Integer>, tsep_src_snk_3pt
				  );
  MeasArgs() { 
    action = ActionType::DWF;
    Ls = 12;
    cfg_stub = "ckpoint_lat";
    rng_stub = "ckpoint_rng";
    gfix_alpha = 0.05;
    gfix_stop_cnd = 1e-9;
    gfix_fourier_accelerate = false;    
    mass = 0.01;
    mobius_scale = 2.0;
    cg_args_exact.algorithm = MixedCGalgorithm::RestartedCG;
    cg_args_exact.tolerance = 1e-8;
    cg_args_exact.restartedcg_inner_tol = 1e-6;
    cg_args_sloppy.algorithm = MixedCGalgorithm::RestartedCG;
    cg_args_sloppy.tolerance = 1e-4;
    cg_args_sloppy.restartedcg_inner_tol = 1e-4;
    nsloppy = 32;
    nexact = 1;
    evec_prec = 2; //2=double 1=single
    tsep_src_snk_3pt = std::vector<Integer>(1, 8);
  }
};

struct Opts{
  std::string arg_file;
  int cfg_start;
  int cfg_lessthan;
  int cfg_step;

  EvecContainerOpts evec_opts_l;
 
  bool load_gfix_cfg = false;
  std::string load_gfix_cfg_stub;

  bool save_gfix_cfg = false;
  std::string save_gfix_cfg_stub;

  bool cps_cfg = false;

  bool unit_gauge = false;
  bool use_evecs = true;

  bool use_split_grid = false;
  std::vector<int> split_grid_proc_layout;

  bool use_gauge_fixing = true;

  Opts(){
    load_gfix_cfg = false;
    save_gfix_cfg = false;
    cps_cfg = false;
    unit_gauge = false;
    use_evecs = true;

    use_split_grid = false;
    use_gauge_fixing = true;
  }
};

std::vector<RealD> operator-(const std::vector<RealD> &a, const std::vector<RealD> &b){
  assert(a.size() == b.size());
  std::vector<RealD> out(a);
  for(int i=0;i<a.size();i++) out[i] -= b[i];
  return out;
}

struct Correlator{
  std::vector<std::vector<RealD> > sloppy; //[tsrc][tsep]
  std::vector<std::vector<RealD> > exact;
  std::vector<RealD> multiplicity; //record multiplicity of exact timeslices for when we source avg
  std::string stub;
  int traj;

  Correlator() = default;
  Correlator(int Lt, const std::string &file_stub, int traj): sloppy(Lt, std::vector<RealD>(Lt,0.)),
							      exact(Lt, std::vector<RealD>(Lt,0.)),
							      stub(file_stub), multiplicity(Lt,0.), traj(traj){}
  void setSloppy(const std::vector<RealD> &cor, int tsrc){
    sloppy[tsrc] = cor;
  }
  void setExact(const std::vector<RealD> &cor, int tsrc, RealD _multiplicity){
    exact[tsrc] = cor;
    multiplicity[tsrc] = _multiplicity;
  }

  void write(){
    asciiWriteArray(sloppy,stub+"_sloppy",traj);
    asciiWriteArray(exact,stub+"_exact",traj);
    asciiWriteArray(multiplicity,stub+"_multiplicity",traj);
  }
  
};


template<typename EvecContainerType>
void invertWithMidProp(LatticePropagatorD &into, LatticePropagatorD &into_mid, const LatticePropagatorD &src, EvecContainerType &eval, ActionsPeriodic &actions, ActionsPeriodic &actions_sub, const MixedCGargs &cg_args, bool use_split_grid){
  if(use_split_grid)
    eval.splitGridMixedPrecInvertWithMidProp(into, into_mid, src, *actions.action_d, *actions_sub.action_d, *actions_sub.action_f, cg_args);
  else
    eval.mixedPrecInvertWithMidProp(into, into_mid, src, *actions.action_d, *actions.action_f, cg_args);
}
template<typename EvecContainerType>
void invert(LatticePropagatorD &into, const LatticePropagatorD &src, EvecContainerType &eval, ActionsPeriodic &actions, ActionsPeriodic &actions_sub, const MixedCGargs &cg_args, bool use_split_grid){
  LatticePropagatorD dummy(src.Grid());
  invertWithMidProp(into,dummy,src,eval,actions,actions_sub,cg_args,use_split_grid);
}
template<typename EvecContainerType>
void invert5D(LatticePropagatorD &into, const LatticePropagatorD &src, EvecContainerType &eval, ActionsPeriodic &actions, ActionsPeriodic &actions_sub, const MixedCGargs &cg_args, bool use_split_grid){
  if(use_split_grid)
    eval.splitGridMixedPrecInvert5D(into, src, *actions.action_d, *actions_sub.action_d, *actions_sub.action_f, cg_args);
  else
    eval.mixedPrecInvert5D(into, src, *actions.action_d, *actions.action_f, cg_args);  
}

template<typename EvecContainerType>
std::vector<RealD> computePionLocalVector3pt(const LatticePropagatorD &prop, int tsrc, int tsep_src_snk, 
			       ActionsPeriodic &actions, ActionsPeriodic &actions_sub,
			       EvecContainerType &evecs, int Lt,  const MixedCGargs &cg_args, bool use_split_grid){
  GridBase* UGridD = prop.Grid();
  
  Gamma gamma4(Gamma::Algebra::GammaT);
  Gamma gamma5(Gamma::Algebra::Gamma5);

  LatticePropagatorD zero_prop(UGridD); zero_prop= Zero();

  int tsnk = (tsrc + tsep_src_snk) % Lt;
  LatticeInteger lcoor(UGridD); LatticeCoordinate(lcoor,3);
  
  LatticePropagatorD seq_src = where( lcoor == tsnk, prop, zero_prop);
  seq_src = gamma5 * seq_src * gamma5; //g5 G(x_snk, tsnk; x_src, tsrc) g5
  
  LatticePropagatorD R_P_seq(UGridD);
  invert(R_P_seq, seq_src, evecs, actions, actions_sub, cg_args, use_split_grid);

  LatticePropagatorD sqb = gamma5 *( adj(prop) * (gamma5 * (gamma4  * R_P_seq))); // [\sum_{x_src'} G(x_src', tsrc;  x,t)] g4 [ \sum_{x_snk,x_src} G(x,t; x_snk, tsnk) g5 G(x_snk, snk; x_src, tsrc) g5 ]

  std::vector<SpinColourMatrixD> Ctm;
  sliceSum(sqb, Ctm, 3);
  assert(Ctm.size() == Lt);
  
  std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;
   
  std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
  for(int t=0;t<Lt;t++)
    Ctx[t] = real( trace(Ctm[t]) );

  std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
  for(int tx=0;tx<Lt;tx++){
    int t = ( tx - tsrc + Lt ) % Lt;
    Ct[t] = Ctx[tx];
  }
  return Ct;
}


template<typename LanczosAction>
void run(const MeasArgs &args, const Opts &opts){
  typedef CayleyFermion5D<WilsonImplD> DefaultAction;
  typedef PeriodicGimplD Gimpl;

  printMem("Program body start");
  
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];

  std::cout << GridLogMessage << "Creating 4d and 5d Grids with Ls=" << args.Ls << std::endl;
  Grids GridsD = makeDoublePrecGrids(args.Ls, latt);
  Grids GridsF = makeSinglePrecGrids(args.Ls, latt);

  printMem("Post Grid creation");
  
  WilsonImplD::ImplParams Params = setupActionParamsPeriodic();

  LatticeGaugeFieldD U_d(GridsD.UGrid);
  LatticeGaugeFieldF U_f(GridsF.UGrid);
  typename LanczosAction::GaugeField & U_lanczos = getInstance<typename LanczosAction::GaugeField>(U_d,U_f);

  printMem("Post gauge field creation");
  
  //Get the actions
  ActionsPeriodic actions(args.action, Params, args.mass, args.mobius_scale, U_d, GridsD, U_f, GridsF);
  printMem("Post action creation");
  
  //Prepare split grid if in use
  std::unique_ptr<Grids> SubGridsD_p, SubGridsF_p;
  ActionsPeriodic actions_sub;
  std::unique_ptr<LatticeGaugeFieldD> U_sub_d_p;
  std::unique_ptr<LatticeGaugeFieldF> U_sub_f_p = nullptr;
  if(opts.use_split_grid){
    printMem("Pre split Grid creation");
    Coordinate proc_sub(opts.split_grid_proc_layout);
    SubGridsD_p.reset(new Grids(makeSplitGrids(GridsD, proc_sub)));
    SubGridsF_p.reset(new Grids(makeSplitGrids(GridsF, proc_sub)));
    printMem("Post split Grid creation");
    U_sub_d_p.reset(new LatticeGaugeFieldD(SubGridsD_p->UGrid));
    U_sub_f_p.reset(new LatticeGaugeFieldF(SubGridsF_p->UGrid));
    printMem("Post split Grid gauge field creation");
    actions_sub = ActionsPeriodic(args.action, Params, args.mass, args.mobius_scale, *U_sub_d_p, *SubGridsD_p, *U_sub_f_p, *SubGridsF_p);
    printMem("Post split-grid Grid/action creation");
  }
  
  //Wall source timeslices
  std::vector<int> src_t(args.nsloppy);
  if(args.nsloppy != 0){ //allow 0 as we might want to just generate and save the eigenvectors / gfix config  
    assert(Lt % args.nsloppy == 0);
    int tsep = Lt / args.nsloppy;
  
    for(int i=0;i<args.nsloppy;i++){    
      int t = i * tsep;
      src_t[i] = t;
    }
    
    assert(args.nexact <= args.nsloppy);
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
    
    //Figure out which timeslices among the sloppy timeslices we will do the AMA correction on
    std::map<int,int> exact_src_t; //timeslice, count   :  we allow for the same slice to be drawn more than once (independent random numbers!)
    for(int i=0;i<args.nexact;i++){
      RealD rt;  random(sRNG,rt); //0->1
      Integer sloppy_idx = (Integer)( rt * args.nsloppy );
      Integer t = src_t[sloppy_idx];

      auto it = exact_src_t.find(t);
      if(it == exact_src_t.end())
	exact_src_t[t] = 1;
      else
	++exact_src_t[t];
    }

    std::cout << GridLogMessage << "Computing with sloppy solves on the following timeslices. The multiplicity of exact solves is given in parentheses." << std::endl;
    int etotal=0;
    for(int s=0;s<args.nsloppy;s++){
      auto eit = exact_src_t.find(src_t[s]);
      int emult = eit == exact_src_t.end() ? 0 : eit->second;
      
      std::cout << src_t[s] << "(" << emult << ") ";
      etotal += emult;
    }
    std::cout << std::endl;
    assert(etotal == args.nexact);

    ///////////////////////////////// Read config /////////////////////////////////////////
    if(opts.unit_gauge){
      std::cout << GridLogMessage << "Setting gauge field to unit gauge" << std::endl;
      SU<Nc>::ColdConfiguration(U_d);
    }else{
      std::cout << GridLogMessage << "Reading configuration" << std::endl;
      opts.cps_cfg ? 
	readCPSconfiguration<Gimpl>(U_d, sRNG, pRNG, traj, args.cfg_stub) :
	readConfiguration<Gimpl>(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);
    }
    printMem("Post configuration read");
    
    ///////////////////////////////// Gauge fix ///////////////////////////////////////////
    if(opts.use_gauge_fixing){
      if(opts.load_gfix_cfg)
	readGaugeFixedConfiguration<Gimpl>(U_d, opts.load_gfix_cfg_stub, traj);    
      else
	CoulombGaugeFix<Gimpl>(U_d, args.gfix_alpha, args.gfix_stop_cnd, args.gfix_fourier_accelerate);
      
      if(opts.save_gfix_cfg)
	writeGaugeFixedConfiguration<Gimpl>(U_d, opts.save_gfix_cfg_stub, traj);
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
      eval.generate(args.lanc_args, traj, *actions.getAction<LanczosAction>(), U_lanczos, pRNG_lanczos, opts.evec_opts_l);
    }else{
      std::cout << GridLogMessage << "Not using light eigenvectors" << std::endl;
    }
    printMem("Post light evec generation");
   
    //////////////////////////////////// Setup outputs //////////////////////////////////////
#define COR(NM) Correlator NM(Lt, #NM, traj)
    COR(C_PP);
    COR(C_J5qP);
    COR(C_PP_WW);
    COR(C_lAP);
    COR(C_cAP);
    
    std::vector<Correlator> C_PVP(args.tsep_src_snk_3pt.size());
    for(int tsepidx=0;tsepidx< C_PVP.size();tsepidx++){
      int tsep = args.tsep_src_snk_3pt[tsepidx];
      std::ostringstream stub; stub << "C_PVP_tsep" << tsep;
      C_PVP[tsepidx] = Correlator(Lt,stub.str(),traj);
    }

    /////////////////////////////////// Trajectory loop //////////////////////////////////////////
    for(int s=0;s<args.nsloppy;s++){ 
      int t0 = src_t[s];
      std::cout << GridLogMessage << "Starting calculation with source timeslice t0=" << t0 << std::endl;
      printMem("Start of timeslice sloppy");
      resetDeviceStackMemory();
      printMem("Post stack memory reset");
      
      //Gauge fixed wall source
      LatticePropagatorD src = wallSource<DefaultAction>(t0, GridsD.UGrid);

      printMem("Post create wall source");

      std::cout << GridLogMessage << "Starting sloppy quark inverse" << std::endl;
      LatticePropagatorD R5D(GridsD.FGrid);
      invert5D(R5D, src, eval, actions, actions_sub, args.cg_args_sloppy, opts.use_split_grid);

      LatticePropagatorD R(GridsD.UGrid), R_mid(GridsD.UGrid);
      convert5Dto4DpropWithMidProp(R,R_mid,R5D, *actions.action_d);

      std::cout << GridLogMessage << "Starting sloppy contractions" << std::endl;
      C_PP.setSloppy( pionCorrelator(t0,R), t0);
      C_J5qP.setSloppy( pionCorrelator(t0,R_mid), t0);
      C_PP_WW.setSloppy( wallSinkPionCorrelator(t0,R), t0);
      C_lAP.setSloppy( pionToLocalAxialTcorrelator(t0,R), t0);
      C_cAP.setSloppy( pionToConservedAxialTcorrelator(t0,R5D,*actions.action_d), t0);
      
      for(int tsepidx=0;tsepidx< C_PVP.size();tsepidx++){
	int tsep = args.tsep_src_snk_3pt[tsepidx];
	std::cout << GridLogMessage << "Starting sloppy PVP 3pt tsep=" << tsep << std::endl;
	C_PVP[tsepidx].setSloppy( computePionLocalVector3pt(R,t0, tsep, actions,actions_sub, eval, Lt, args.cg_args_sloppy, opts.use_split_grid),   t0 );
      }
      
      auto eit = exact_src_t.find(t0);
      if(eit != exact_src_t.end()){
	int multiplicity = eit->second;

	//Do exact solve
	std::cout << GridLogMessage << "Starting exact quark inverse" << std::endl;
	invert5D(R5D, src, eval, actions, actions_sub, args.cg_args_exact, opts.use_split_grid);
	convert5Dto4DpropWithMidProp(R,R_mid,R5D, *actions.action_d);

	std::cout << GridLogMessage << "Starting exact contractions" << std::endl;
	C_PP.setExact( pionCorrelator(t0,R), t0, multiplicity);
	C_J5qP.setExact( pionCorrelator(t0,R_mid), t0, multiplicity);
	C_PP_WW.setExact( wallSinkPionCorrelator(t0,R), t0, multiplicity);
	C_lAP.setExact( pionToLocalAxialTcorrelator(t0,R), t0, multiplicity);
	C_cAP.setExact( pionToConservedAxialTcorrelator(t0,R5D,*actions.action_d), t0, multiplicity);
	
	for(int tsepidx=0;tsepidx< C_PVP.size();tsepidx++){
	  int tsep = args.tsep_src_snk_3pt[tsepidx];
	  std::cout << GridLogMessage << "Starting exact PVP 3pt tsep=" << tsep << std::endl;
	  C_PVP[tsepidx].setExact( computePionLocalVector3pt(R,t0, tsep, actions,actions_sub, eval, Lt, args.cg_args_exact, opts.use_split_grid),   t0, multiplicity );
	}
      }
      
      C_PP.write();
      C_J5qP.write();
      C_PP_WW.write();
      C_lAP.write();
      C_cAP.write();
      for(int tsepidx=0;tsepidx< C_PVP.size();tsepidx++)
	C_PVP[tsepidx].write();
    }
    
    
    
  }//traj
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
    }else if(sargv == "--load_gfix_cfg"){
      opts.load_gfix_cfg=true;
      opts.load_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--save_gfix_cfg"){
      opts.save_gfix_cfg=true;
      opts.save_gfix_cfg_stub = argv[i+1];
    }else if(sargv == "--disable_evecs"){
      opts.use_evecs = false;
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
    run<CayleyFermion5D<WilsonImplD> >(args,opts);   
  }else if(args.evec_prec == 1){
    run<CayleyFermion5D<WilsonImplF> >(args,opts);
  }else{
    assert(0);
  }
   
  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
