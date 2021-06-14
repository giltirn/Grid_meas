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
				  double, ml,
				  double, cg_stop,
				  double, cg_stop_inner,
				  int, nsrc,
				  LanczosParameters, lanc_args
				  );
  MeasArgs() { 
    action = ActionType::DWF;
    Ls = 12;
    GparityDirs = {1,0,0};
    cfg_stub = "ckpoint_lat";
    rng_stub = "ckpoint_rng";
    gfix_alpha = 0.05;
    ml = 0.01;
    mobius_scale = 2.0;
    cg_stop = 1e-8;
    cg_stop_inner = 1e-6;
    nsrc = 32;
  }
};

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  assert(argc >= 5);
  std::string arg_file = argv[1];
  int cfg_start = std::stoi(argv[2]);
  int cfg_lessthan = std::stoi(argv[3]);
  int cfg_step = std::stoi(argv[4]);
 
  bool load_evecs = false;
  std::string load_evecs_stub, load_evals_stub;

  bool save_evecs = false;
  std::string save_evecs_stub, save_evals_stub;

  bool cps_cfg = false;
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

  auto UGridD   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  auto UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(args.Ls, UGridD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(args.Ls, UGridD);

  auto UGridF   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  auto UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(args.Ls, UGridF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(args.Ls, UGridF);

  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = args.GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = args.GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction


  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(UGridD); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD U_d(UGridD);
  LatticeGaugeFieldF U_f(UGridF);

  CayleyFermion5D<GparityWilsonImplD>* action_d = createActionD(args.action, Params, args.ml, args.mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<GparityWilsonImplF>* action_f = createActionF(args.action, Params, args.ml, args.mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  //Not going to worry about a rotationally invariant operator,
  //just use p_j = pi/2L for both quark and antiquark in direction j
  std::vector<int> p1(3,0);
  for(int i=0;i<3;i++)
    p1[i] = args.GparityDirs[i] ? 1 : 0;
  std::vector<double> p1_phys = getPhysicalMomentum(p1);
  
  std::vector<int> p2 = p1;
  std::vector<double> p2_phys = p1_phys;

  //Z2wall source
  assert(Lt % args.nsrc == 0);
  int tsep = Lt / args.nsrc;
  std::vector<int> src_t(args.nsrc);
  
  std::vector<LatticeSCFmatrixD> sources(args.nsrc, UGridD);
  for(int i=0;i<args.nsrc;i++){    
    int t = i * tsep;
    src_t[i] = t;
    sources[i] = Z2wallSource(t, pRNG, UGridD);
  }

  LatticeComplexD p1_src_phase_field = phaseField(p1_phys, UGridD); //exp(-i \vec p \cdot \vec x)
  p1_src_phase_field = conjugate(p1_src_phase_field); //exp(+i \vec p \cdot \vec x)

  const LatticeComplexD &p2_src_phase_field = p1_src_phase_field;

  std::vector<RealD> eval;
  std::vector<FermionFieldD> evec;

  //Start traj loop
  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    cps_cfg ? 
      readCPSconfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub) :
      readConfiguration(U_d, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);

    precisionChange(U_f, U_d);

    action_d->ImportGauge(U_d);
    action_f->ImportGauge(U_f);

    if(load_evecs){
      readEigenvalues(eval, evec, FrbGridD, load_evals_stub, load_evecs_stub, traj);
    }else{
      computeEigenvalues<CayleyFermion5D<GparityWilsonImplD>, FermionFieldD>(eval, evec, args.lanc_args, FGridD, FrbGridD, U_d, *action_d, pRNG);
    }
    if(save_evecs) saveEigenvalues(eval, evec, save_evals_stub, save_evecs_stub, traj);

    std::vector<RealD> Ct_pion(Lt, 0);
    std::vector<RealD> Ct_j5q(Lt, 0);

    for(int s=0;s<args.nsrc;s++){ 
      int t0 = src_t[s];
      std::cout << GridLogMessage << "Starting calculation with source timeslice t0=" << t0 << std::endl;

      LatticeSCFmatrixD src_p1 = p1_src_phase_field * sources[s];

      LatticeSCFmatrixD Rp1(UGridD), Rp1_mid(UGridD);
      mixedPrecInvertWithMidProp(Rp1, Rp1_mid, src_p1, *action_d, *action_f, args.cg_stop, args.cg_stop_inner, &eval, &evec);
      
      const LatticeSCFmatrixD &Rp2 = Rp1;
      const LatticeSCFmatrixD &Rp2_mid = Rp1_mid;

      std::vector<RealD> Cts_pion = momWallSourcePionCorrelator(p1, p2, t0, Rp1, Rp2);
      std::vector<RealD> Cts_j5q = momWallSourcePionCorrelator(p1, p2, t0, Rp1_mid, Rp2_mid);

      for(int t=0;t<Lt;t++){ //average over sources
	Ct_pion[t] += Cts_pion[t] / RealD(args.nsrc);
	Ct_j5q[t] += Cts_j5q[t] / RealD(args.nsrc);
      }

    }

    asciiWriteArray(Ct_pion, "pion_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
    asciiWriteArray(Ct_j5q, "j5q_mom" + momstr(p1) + "_mom" + momstr(p2), traj);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
