#include<Grid/Grid.h>

using namespace Grid;

typedef GparityWilsonImplD::FermionField FermionField;
typedef GparityWilsonImplD::SiteSpinor SiteSpinor;

GRID_SERIALIZABLE_ENUM(ActionType, undef, DWF, 1, Mobius, 2);

struct MeasArgs: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MeasArgs,
				  ActionType, action,
				  int, Ls,
				  double, mobius_scale,
				  std::vector<Integer>, GparityDirs,
				  std::string, cfg_stub,
				  std::string, rng_stub,
				  double, gfix_alpha,
				  double, ml
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
  }
};

bool fileExists(const std::string &fn){
  std::ifstream f(fn);
  return f.good();
}


void readConfiguration(LatticeGaugeFieldD &U,
		       GridSerialRNG &sRNG,
		       GridParallelRNG &pRNG,
		       int traj, 
		       const std::string &cfg_stub,
		       const std::string &rng_stub){

  CheckpointerParameters p(cfg_stub, rng_stub);
  NerscHmcCheckpointer<ConjugateGimplD> cp(p);

  cp.CheckpointRestore(traj, U, sRNG, pRNG);
}

//For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
void readCPSconfiguration(LatticeGaugeFieldD &U,
			  GridSerialRNG &sRNG,
			  GridParallelRNG &pRNG,
			  int traj, 
			  const std::string &cfg_stub){

  NerscIO::exitOnReadPlaquetteMismatch() = false;

  CheckpointerParameters p(cfg_stub, "pooh");
  NerscHmcCheckpointer<ConjugateGimplD> cp(p);
  typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
    
  std::string config, dummy;
  cp.build_filenames(traj, p, config, dummy);  
  cp.check_filename(config);
  
  FieldMetaData header;
  NerscIO::readConfiguration<GaugeStats>(U, header, config);

  NerscIO::exitOnReadPlaquetteMismatch() = true;

  std::vector<int> seeds4({traj, traj+1, traj+2, traj+3});
  pRNG.SeedFixedIntegers(seeds4);
  sRNG.SeedFixedIntegers(seeds4); 
}
  


typedef iMatrix<iMatrix<iMatrix<vComplexD, Nc>, Ns>, Ngp> vSCFmatrix;
typedef iMatrix<iMatrix<iMatrix<ComplexD, Nc>, Ns>, Ngp> SCFmatrix;

typedef Lattice<vSCFmatrix> LatticeSCFmatrix;

SCFmatrix unitSiteSrc(){
  SCFmatrix out;
  out = Zero();
  for(int f=0;f<Ngp;f++)
  for(int s=0;s<Ngp;s++)  
  for(int c=0;c<Ngp;c++)  
    out(f,f)(s,s)(c,c) = 1.;
  return out;
}

LatticeSCFmatrix pointSource(const Coordinate &coord, GridBase* UGrid){
  SCFmatrix one = unitSiteSrc();
  
  LatticeSCFmatrix src(UGrid);
  src = Zero();
  pokeSite(one, src, coord);
  return src;
}

//From the matrices extract the column with the given flavor, spin and color index
FermionField extractColumn(const LatticeSCFmatrix &from, const int fs, const int ss, const int cs){
  FermionField into(from.Grid());

  typedef FermionField::vector_type vector_type;
  Lattice<iScalar<iScalar<iVector<vector_type, Nc> > > > CV(from.Grid());
  Lattice<iScalar<iVector<iVector<vector_type, Nc>, Ns> > > SCV(from.Grid());
  
  for(int fv=0;fv<Ngp;fv++){
	  
    auto elem_ff = PeekIndex<GparityFlavourIndex>(from, fv, fs);
	  
    for(int sv=0;sv<Ns;sv++){
	    
      auto elem_ff_ss = PeekIndex<SpinIndex>(elem_ff, sv, ss);
	    
      for(int cv=0;cv<Nc;cv++){    
	      
	auto elem_ff_ss_cc = PeekIndex<ColourIndex>(elem_ff_ss, cv, cs); //scalar<scalar<scalar> > >
	      
	PokeIndex<ColourIndex>(CV, elem_ff_ss_cc, cv);
      }
      PokeIndex<SpinIndex>(SCV, CV, sv);
    }
    PokeIndex<GparityFlavourIndex>(into, SCV, fv);
  }
  return into;
}


void insertColumn(LatticeSCFmatrix &into, const FermionField &from, const int fs, const int ss, const int cs){
  typedef FermionField::vector_type vector_type;
  Lattice<iScalar<iMatrix<iMatrix<vector_type, Nc>, Ns> > > SCM(into.Grid());
  Lattice<iScalar<iScalar<iMatrix<vector_type, Nc> > > > CM(into.Grid());
 
  for(int fv=0;fv<Ngp;fv++){
	  
    auto elem_f = PeekIndex<GparityFlavourIndex>(from, fv);
	  
    for(int sv=0;sv<Ns;sv++){
	    
      auto elem_f_s = PeekIndex<SpinIndex>(elem_f, sv);
	    
      for(int cv=0;cv<Nc;cv++){    
	      
	auto elem_f_s_c = PeekIndex<ColourIndex>(elem_f_s, cv); //scalar<scalar<scalar> > >
	      
	PokeIndex<ColourIndex>(CM, elem_f_s_c, cv, cs);
      }
      PokeIndex<SpinIndex>(SCM, CM, sv, ss);
    }
    PokeIndex<GparityFlavourIndex>(into, SCM, fv, fs);
  }
}




template<typename FermionAction>
LatticeSCFmatrix invert(const LatticeSCFmatrix &msrc, FermionAction &action){
  ConjugateGradient<FermionField> CG(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField> solver(CG);
  
  GridBase* FGrid = action.FermionGrid();

  FermionField src_5d(FGrid);
  FermionField sol_5d(FGrid);
  
  FermionField sol_4d(msrc.Grid());
  LatticeSCFmatrix msol(msrc.Grid());

  //Columns of the matrix are the source vectors
  for(int f=0;f<Ngp;f++){
    for(int s=0;s<Ns;s++){
      for(int c=0;c<Nc;c++){
	FermionField src_4d = extractColumn(msrc, f,s,c);
	action.ImportPhysicalFermionSource(src_4d, src_5d);

	solver(action, src_5d, sol_5d);
	
	action.ExportPhysicalFermionSource(sol_5d, sol_4d);
	
	insertColumn(msol, sol_5d, f,s,c);
      }
    }
  }
  return msol;
}


CayleyFermion5D<GparityWilsonImplD>* createAction(ActionType action,
						  const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						  LatticeGaugeField &Umu,
						  GridCartesian         &FiveDimGrid,
						  GridRedBlackCartesian &FiveDimRedBlackGrid,
						  GridCartesian         &FourDimGrid,
						  GridRedBlackCartesian &FourDimRedBlackGrid
						  ){
  double bpc = mobius_scale;
  double bmc = 1.0;
  double b = (bpc + bmc)/2.;
  double c = (bpc - bmc)/2.;
    
  switch(action){
  case ActionType::DWF:
    return new GparityDomainWallFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, Params);
  case ActionType::Mobius:
    return new GparityMobiusFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, b, c, Params);
  };
  return nullptr;
};

LatticeSCFmatrix invert(const LatticeSCFmatrix &msrc, CayleyFermion5D<GparityWilsonImplD>* action){
  return invert(msrc, *action);
}




//f(\vec x, t) = exp(-i \vec p \cdot \vec x)
LatticeComplexD phaseField(const std::vector<double> &p, GridBase* Grid){
  LatticeComplexD out(Grid);
  LatticeComplexD coor(Grid);
  out=Zero();
  for(int mu=0;mu<3;mu++){
    LatticeCoordinate(coor,mu);
    out = out + p[mu] * coor;
  }
  ComplexD ci(0.0,-1.0);
  out = exp(out*ci);
  return out;
}

//exp(-i \vec p \cdot x0)
ComplexD phase(const std::vector<double> &p, const Coordinate &x0){
  ComplexD out(0);
  for(int mu=0;mu<3;mu++){
    out = out + p[mu] * x0[mu];
  }
  ComplexD ci(0.0,-1.0);
  out = exp(out*ci);
  return out;
}

//Get the physical momentum in lattice units. Input units are :  2pi/L  (non-Gparity)   pi/2L (G-parity)
std::vector<double> getPhysicalMomentum(const std::vector<int> &p){
  assert(p.size() == 3);
  std::vector<int> gpdirs = ConjugateGimplD::getDirections();
  Coordinate latt = GridDefaultLatt();
  std::vector<double> pphys(3,0);
  for(int i=0;i<3;i++){
    double punit = gpdirs[i] ? M_PI/2./latt[i] : 2*M_PI/latt[i];
    pphys[i] = p[i]*punit;
  }
  return pphys;
}


//C(t) =   \sum_{\vec x, \vec y_1, \vec y_2} e^{-i( \vec p\cdot \vec x - \vec p_1 \cdot \vec y_1 - \vec p_2 \cdot \vec y_2 ) } tr [ \sigma_3 G(\vec x, \vec y_1) \sigma_3 \phi(\vec y_1, \vec y_2)  G(\vec y_2, \vec x) ]
//\vec p_1 + \vec p_2 = \vec p

//For point sources:
//C(t) =   \sum_{\vec x} e^{-i\vec p\cdot( \vec x - \vec y_0 ) } tr [ \sigma_3 G(\vec x, \vec y_0) \sigma_3 G(\vec y_0, \vec x) ]
//     =   \sum_{\vec x} e^{-i\vec p\cdot( \vec x - \vec y_0 ) } tr [ \sigma_3 G(\vec x, \vec y_0) \sigma_3 \gamma^5 G^\dagger(\vec x, \vec y_0) \gamma^5 ]
//     =   \sum_{\vec x} e^{-i\vec p\cdot \vec x}e^{i\vec p\cdot \vec y_0} tr [ \sigma_3 G(\vec x, \vec y_0) \sigma_3 \gamma^5 G^\dagger(\vec x, \vec y_0) \gamma^5 ]

//No need for source flavor projection as not a smeared source

std::vector<RealD> pointSourcePionCorrelator(const std::vector<int> &p, const Coordinate &y0, const LatticeSCFmatrix& G){
  std::vector<double> pphys = getPhysicalMomentum(p);
  ComplexD src_phase = conj( phase(pphys,y0 ) ); //exp(i \vec p \cdot y_0 )
  LatticeComplexD snk_phase_field = phaseField(pphys, G.Grid()); //exp(-i \vec p \cdot \vec x)
  
  GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
  Gamma gamma5(Gamma::Algebra::Gamma5);

  LatticeSCFmatrix sqb =  snk_phase_field *(sigma3* (gamma5 * G * sigma3 * (gamma5 *  adj(G)))); //brutally force the compiler to pair up the operands in a way that makes sense
  
  //LatticeSCFmatrix sqb = snk_phase_field * sigma3 * gamma5 * G * sigma3 * gamma5 * adj(G);
  std::vector<SCFmatrix> Ctm;
  sliceSum(sqb, Ctm, 3);
  
  int Lt = Ctm.size();
  std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
  for(int t=0;t<Lt;t++)
    Ctx[t] = real( trace(Ctm[t]) * src_phase );

  std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
  for(int tx=0;tx<Lt;tx++){
    int t = ( tx - y0[3] + Lt ) % Lt;
    Ct[t] = Ctx[tx];
  }
  return Ct;
}


void write(const std::vector<RealD> &data, const std::string &stub, const int traj){
  std::string filename = stub + "." + std::to_string(traj);
  std::ofstream of(filename);

  for(int t=0;t<data.size();t++){
    of << data[t] << (t < data.size() - 1 ? " " : "");
  }
  of << std::endl;
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

  auto UGridD   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  auto UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(args.Ls, UGridD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(args.Ls, UGridD);

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

  LatticeGaugeField U(UGridD);

  CayleyFermion5D<GparityWilsonImplD>* action = createAction(args.action, Params, args.ml, args.mobius_scale, U, *FGridD, *FrbGridD, *UGridD, *UrbGridD);


  std::vector<LatticeInteger> coor(4, UGridD);
  for(int i=0;i<4;i++) LatticeCoordinate(coor[i], i);

  Coordinate pnt({0,0,0,0});
  LatticeSCFmatrix point_src = pointSource(pnt, UGridD);

  std::vector<int> p(3,0);
  for(int i=0;i<3;i++)
    p[i] = args.GparityDirs[i] ? 2 : 0;


  for(int traj = cfg_start; traj < cfg_lessthan; traj += cfg_step){
    cps_cfg ? 
      readCPSconfiguration(U, sRNG, pRNG, traj, args.cfg_stub) :
      readConfiguration(U, sRNG, pRNG, traj, args.cfg_stub, args.rng_stub);

    action->ImportGauge(U);

    LatticeSCFmatrix G = invert(point_src, action);

    std::vector<RealD> Ct = pointSourcePionCorrelator(p, pnt, G);

    write(Ct, "pion_111", traj);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
