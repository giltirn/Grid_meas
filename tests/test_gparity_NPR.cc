#include<common/propagator_invert.h>
#include<common/sources.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/momentum.h>

using namespace Grid;
using namespace GridMeas;

namespace NPR{
  typedef std::array<int,4> Momentum; //units of pi/2L in x,y,z,  pi/L in t  (to accommodate aprd/prd)

  //n_pj = L pj /pi - 1/2   for pj in lattice units
  //for  Pj in units of pi/2L
  //n_pj = L/pi ( Pj pi / 2L ) - 1/2 = Pj/2 - 1/2
  int _2npj(const Momentum &p, const int j){ return p[j] - 1; }

  //Can be negative!
  bool isMultipleOfFour(int i){
    return (i + 400) % 4 == 0;
  }

  //NOTE: for BCs vector of int, the first 3 directions indicate G-parity dirs and the fourth if the time direction is APRD. Values are restricted to 0 (off), 1 (on)
  
  //Allowed momenta differ by 2pi n/L = 4n (pi/2L)   in G-parity directions for integer n
  bool allowedDirection(const Momentum &p, const std::vector<int> &BCs){
    if(BCs[0] && BCs[1] && BCs[2]){ 
      return isMultipleOfFour(p[1]-p[0]) && isMultipleOfFour(p[2]-p[0]) && isMultipleOfFour(p[2]-p[1]);
    }else if(BCs[0] && BCs[1]){
      return isMultipleOfFour(p[1]-p[0]);
    }else if(BCs[0] && BCs[2]){
      return isMultipleOfFour(p[2]-p[0]);
    }else if(BCs[1] && BCs[2]){
      return isMultipleOfFour(p[2]-p[1]);
    }else{
      return true; //GP in 1 or 0 dirs
    }
  }

  Momentum operator-(const Momentum &p1, const Momentum &p2){
    Momentum out;
    for(int i=0;i<4;i++) out[i] = p1[i] - p2[i];
    return out;
  }
  
  std::vector<double> momUnits(int L, int Lt, const std::vector<int> &BCs){
    std::vector<double> units(4);
    for(int i=0;i<3;i++) units[i] = BCs[i] ? M_PI/(2*L) : 2*M_PI/L;
    units[3] = BCs[3] ? M_PI/Lt : 2*M_PI/Lt;
    return units;
  }

  //Return the squared magnitude in lattice units
  double sqMag(const Momentum &p, const std::vector<double> &units){
    double out = 0.;
    for(int i=0;i<4;i++)
      out += pow(p[i]*units[i], 2);
    return out;
  }

  std::vector<double> toPhase(const Momentum &p, const std::vector<double> &units){
    std::vector<double> out(4);
    for(int i=0;i<4;i++) out[i] = p[i] * units[i];
    return out;
  }

  std::ostream & operator<<(std::ostream &os, const Momentum &p){
    os << "(" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ")";
    return os;
  }

  //Get the first allowed symmetric momentum combination with mu>mumin, returns mu
  double getSymmetricMomentumComb(std::vector<double> &p1_out, std::vector<double> &p2_out, const double ainv, const double mumin, int L, int Lt, const std::vector<int> &BCs, const int mmax=7){
    std::vector<double> units = momUnits(L,Lt,BCs);

    std::cout << "Getting allow momenta" << std::endl;
    std::vector<Momentum> allowed;

    Momentum p = {0,0,0,0};
    for(int m0=-mmax;m0<=mmax;m0++){
      p[0] = BCs[0] ? 2*m0 + 1 : m0; //units of pi/2L or 2pi/L
      for(int m1=-mmax;m1<=mmax;m1++){
	p[1] = BCs[1] ? 2*m1 + 1 : m1;
	for(int m2=-mmax;m2<=mmax;m2++){
	  p[2] = BCs[2] ? 2*m2 + 1 : m2;
	  for(int m3=-mmax;m3<=mmax;m3++){
	    p[3] = BCs[3] ? 2*m3 + 1 : m3; //units of pi/L or 2pi/L

	    if(allowedDirection(p,BCs)){
	      allowed.push_back(p);
	    }
	  }
	}
      }
    }
    std::cout << "Searching for allowed momentum combinations using " << allowed.size() << " allowed momenta" << std::endl;
  
    double best_mu = 1e20;
    bool found_one = false;
    int n_symm = 0;
    double mu_lo = 1e20;
    double mu_hi = 0;

    for(auto const &p1 : allowed){
      for(auto const &p2 : allowed){    
	double p1sq = sqMag(p1,units);
	double p2sq = sqMag(p2,units);
	double p1mp2sq = sqMag(p1-p2,units);
      
	if( fabs(p1sq-p2sq)<1e-5 && fabs(p1mp2sq - p2sq)<1e-5 ){
	  double mu = sqrt(p1sq)*ainv;
	  ++n_symm;
	  mu_lo = std::min(mu_lo, mu);
	  mu_hi = std::max(mu_hi, mu);

	  if(mu >= mumin && mu < best_mu){ //try to get as close as possible
	    std::cout << "p1:" << p1 << "  p2:" << p2 << " (p1-p2):" << (p1-p2) << " |p1|^2:" << p1sq << " |p2|^2:" << p2sq << " |p1-p2|^2:" << p1mp2sq << " mu: " << mu << " GeV" << std::endl;
	    p1_out = toPhase(p1, units);
	    p2_out = toPhase(p2, units);
	    best_mu = mu;
	    found_one = true;
	  }
	}
      }
    }
    if(n_symm > 0) std::cout << "Completed search of " << allowed.size()*allowed.size() << " combinations, found " << n_symm << " symmetric combinations with mu in range " << mu_lo << "->" << mu_hi << " GeV" << std::endl;
    else std::cout << "Completed search of " << allowed.size()*allowed.size() << " combinations, found NO symmetric combinations" << std::endl;

    if(found_one){
      std::cout << "Found a symmetric combination with mu=" << best_mu << " GeV, close to target of " << mumin << " GeV" << std::endl;
      return best_mu;
    }else{
      std::cout << "Found NO symmetric combination with mu>=" << mumin << " GeV" << std::endl;
      return -1;
    }
  }
};
  
inline std::vector<double> minusMom(const std::vector<double> &p){ return std::vector<double>({-p[0],-p[1],-p[2],-p[3]}); }
  
template<typename Field>
typename Field::scalar_object sinkTransform(const std::vector<double> &p, const Field &field){
  LatticeComplexD pfield = phaseFieldFour(p, field.Grid());
  Field tmp = pfield * field;
  return sum(tmp);
}

typedef typename LatticeSCFmatrixD::scalar_object SCFmatrixD;

SCFmatrixD invert(const SCFmatrixD &m){
  Eigen::MatrixXcd mm = Eigen::MatrixXcd::Zero(24,24);
  for(int f1=0;f1<2;f1++){
    for(int s1=0;s1<4;s1++){
      for(int c1=0;c1<3;c1++){
	int scf1 = c1+3*(s1 + 4*f1);
	for(int f2=0;f2<2;f2++){
	  for(int s2=0;s2<4;s2++){
	    for(int c2=0;c2<3;c2++){	      
	      int scf2 = c2+3*(s2 + 4*f2);

	      mm(scf1,scf2) = m(f1,f2)(s1,s2)(c1,c2);
	    }
	  }
	}
      }
    }
  }
  Eigen::MatrixXcd mminv = mm.inverse();
  
  SCFmatrixD out;
  for(int f1=0;f1<2;f1++){
    for(int s1=0;s1<4;s1++){
      for(int c1=0;c1<3;c1++){
	int scf1 = c1+3*(s1 + 4*f1);
	for(int f2=0;f2<2;f2++){
	  for(int s2=0;s2<4;s2++){
	    for(int c2=0;c2<3;c2++){	      
	      int scf2 = c2+3*(s2 + 4*f2);

	      out(f1,f2)(s1,s2)(c1,c2) = mminv(scf1,scf2);
	    }
	  }
	}
      }
    }
  }
  return out;
}

template<typename Field>
Field mulF11left(const Field &f){
  static GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
  Field out = ComplexD(0.5)*( f + sigma3 * f );
  return out;
}
template<typename Field>
Field mulF11right(const Field &f){
  static GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
  Field out = ComplexD(0.5)*( f + f* sigma3 );
  return out;
}
template<typename Field>
Field mulF11leftright(const Field &f){
  static GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
  Field out = ComplexD(0.25)*( f + sigma3 * f + f * sigma3 + sigma3 * (f * sigma3) );
  return out;
}



SpinColourMatrixD invert(const SpinColourMatrixD &m){
  Eigen::MatrixXcd mm = Eigen::MatrixXcd::Zero(12,12);
  for(int s1=0;s1<4;s1++){
    for(int c1=0;c1<3;c1++){
      int sc1 = c1+3*s1;
      for(int s2=0;s2<4;s2++){
	for(int c2=0;c2<3;c2++){	      
	  int sc2 = c2+3*s2;

	  mm(sc1,sc2) = m()(s1,s2)(c1,c2);
	}
      }
    }
  }
  Eigen::MatrixXcd mminv = mm.inverse();
  
  SpinColourMatrixD out;
  for(int s1=0;s1<4;s1++){
    for(int c1=0;c1<3;c1++){
      int sc1 = c1+3*s1;
      for(int s2=0;s2<4;s2++){
	for(int c2=0;c2<3;c2++){	      
	  int sc2 = c2+3*s2;

	  out()(s1,s2)(c1,c2) = mminv(sc1,sc2);
	}
      }
    }
  }
  return out;
}


int get_n_GP(const std::vector<double> &p, const int L, const std::vector<int> &BCs){
  int first_GPdir = -1;
  for(int i=0;i<3;i++) if(BCs[i]){ first_GPdir = i; break; }
  if(first_GPdir == -1) assert(0 && "get_n_GP: Using twisted BCs for the comparison requires at least 1 G-parity direction");

  int n = int( floor(  (p[first_GPdir] * 2*L/M_PI - 1.)/2. + 0.5 ) );
  if( fabs( (2*n+1)*M_PI/2./L - p[first_GPdir] ) > 1e-05 ) assert(0 && "get_n_GP: Failed calc check");
  return n;
}


//Checks:
//1) GPBC in 0-dirs vs periodic, random gauge field, identical momenta:  results identical

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  assert(argc >= 4);

  std::vector<int> BCs; //x,y,z: 1 for GP, 0 for prd;   t: 1 for aprd, 0 for prd
  GridCmdOptionIntVector(argv[1], BCs);
  assert(BCs.size() == 4);

  int nGPdirs = 0;
  for(int i=0;i<3;i++) nGPdirs += BCs[i];

  int tBC_P; std::string a2(argv[2]);//grr
  GridCmdOptionInt(a2, tBC_P); //0=prd, 1=aprd

  double mumin; std::string a3(argv[3]);
  GridCmdOptionFloat(a3, mumin);

  int mmax=7;
  double ainv = 1.73;
  int Ls = 12;
  RealD ml = 0.05;
  int Gamma_idx = 0; // 0: scalar, 15: pseudoscalar, 1-3: gamma_mu
 
  for(int i=0;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "-ainv"){
      assert(i<argc-1);
      std::string narg(argv[i+1]);
      GridCmdOptionFloat(narg, ainv);
      std::cout << "Set ainv to " << ainv << std::endl;
    }else if(sarg == "-mmax"){
      assert(i<argc-1);
      std::string narg(argv[i+1]);
      GridCmdOptionInt(narg, mmax);
      std::cout << "Set mmax to " << mmax << std::endl;
    }else if(sarg == "-Ls"){
      assert(i<argc-1);
      std::string narg(argv[i+1]);
      GridCmdOptionInt(narg, Ls);
      std::cout << "Set Ls to " << Ls << std::endl;
    }else if(sarg == "-ml"){
      assert(i<argc-1);
      std::string narg(argv[i+1]);
      GridCmdOptionFloat(narg, ml);
      std::cout << "Set ml to " << ml << std::endl;
    }else if(sarg == "-Gamma"){
      assert(i<argc-1);
      std::string narg(argv[i+1]);
      GridCmdOptionInt(narg, Gamma_idx);
      std::cout << "Set Gamma to " << Gamma_idx << std::endl;
    }
  }

  std::map<int, Gamma::Algebra> gmap = { {0,Gamma::Algebra::Identity}, {15,Gamma::Algebra::Gamma5}, 
					 {1,Gamma::Algebra::GammaX}, {2,Gamma::Algebra::GammaY}, {3,Gamma::Algebra::GammaZ}, {4,Gamma::Algebra::GammaT} };

  assert(gmap.count(Gamma_idx) && "Unknown choice of Gamma");
  Gamma GammaOp(gmap[Gamma_idx]);


  Coordinate latt = GridDefaultLatt();
  std::cout << "Lattice size " << latt << std::endl;

  assert(latt[0] == latt[1] && latt[0] == latt[2]);
  int Lt = latt[3];
  int L = latt[0];
  size_t V3d = L*L*L;
  size_t V4d = V3d*Lt;

  //Get the quark momenta
#define TWISTED_PRD


  std::cout << "mu_min = " << mumin << " GeV, a^-1 = " << ainv << " GeV" << std::endl;
  
  std::vector<double> p1_GP, p2_GP;
  double mu_GP = NPR::getSymmetricMomentumComb(p1_GP, p2_GP, ainv, mumin, L, Lt, BCs, mmax);
  if(mu_GP != -1.) std::cout << "Found symmetric momentum combination for Gparity lattice with mu=" << mu_GP << std::endl;

#ifndef TWISTED_PRD
  std::vector<double> p1_P, p2_P;
  std::vector<int> BCs_P = {0,0,0,tBC_P};
  double mu_P = NPR::getSymmetricMomentumComb(p1_P, p2_P, ainv, mumin, L, Lt, BCs_P, mmax);
  if(mu_P != -1.) std::cout << "Found symmetric momentum combination for periodic lattice with mu=" << mu_P << std::endl;
#else
  double mu_P = mu_GP;
  std::vector<double> p1_P(p1_GP), p2_P(p2_GP);
#endif

  if(mu_GP==-1. || mu_P==-1.) return 0;

  auto UGridD   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  auto UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridD);

  auto UGridF   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  auto UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridF);

  assert(Nd == 4);
  std::vector<int> gauge_BCs({ BCs[0], BCs[1], BCs[2], 0 });
  ConjugateGimplD::setDirections(gauge_BCs); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  Params.twists = Coordinate(BCs);

  std::cout << "Using fermion BCs " << BCs << " and gauge BCs " << gauge_BCs << std::endl;

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(UGridD); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD U_d(UGridD);
  LatticeGaugeFieldF U_f(UGridF);
  
  //Unit gauge
  LorentzColourMatrixD SU3one;
  SU3one = Zero();
  for(int d=0;d<Nd;d++)
    for(int c=0;c<Nc;c++)  
      SU3one(d)()(c,c) = 1.;
  U_d = SU3one;
  
  //gaussian(pRNG, U_d);
  precisionChange(U_f, U_d);
  
  ActionType action = ActionType::Mobius;
  RealD mobius_scale = 2.0;

  CayleyFermion5D<GparityWilsonImplD>* action_d = createActionD(action, Params, ml, mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<GparityWilsonImplF>* action_f = createActionF(action, Params, ml, mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  
#ifdef TWISTED_PRD
  WilsonImplD::ImplParams PParams;
  PParams.boundary_phases[3] = BCs[3] ? -1. : 1.; //same as GP

  int twist_sign_p1=0, twist_sign_p2=0;

  if(nGPdirs > 0){
    //For G-parity momentum  p = (2n+1)\pi / 2L,   we need different twists depending on whether n is odd or even
    int n_p1 = get_n_GP(p1_GP, L, BCs);
    int n_p2 = get_n_GP(p2_GP, L, BCs);
    
    twist_sign_p1 = abs(n_p1) % 2 == 0 ? 1 : -1;
    twist_sign_p2 = abs(n_p2) % 2 == 0 ? 1 : -1;
    
    std::cout << "n_p1=" << n_p1 << " twist sign " << twist_sign_p1 << ",  n_p2=" << n_p2 << " twist sign " << twist_sign_p2 << std::endl;
  }

  for(int mu=0;mu<3;mu++)
    if(BCs[mu])  PParams.boundary_phases[mu] = ComplexD( std::polar(1., twist_sign_p1 * M_PI/2. ) );
 
  CayleyFermion5D<WilsonImplD>* paction_d_p1 = createPeriodicActionD(action, PParams, ml, mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<WilsonImplF>* paction_f_p1 = createPeriodicActionF(action, PParams, ml, mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  for(int mu=0;mu<3;mu++)
    if(BCs[mu])  PParams.boundary_phases[mu] = ComplexD( std::polar(1., twist_sign_p2 * M_PI/2. ) ); 

  CayleyFermion5D<WilsonImplD>* paction_d_p2 = createPeriodicActionD(action, PParams, ml, mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<WilsonImplF>* paction_f_p2 = createPeriodicActionF(action, PParams, ml, mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

#else
  WilsonImplD::ImplParams PParams;
  PParams.boundary_phases[3] = tBC_P ? -1 : 1;

  CayleyFermion5D<WilsonImplD>* paction_d_p1 = createPeriodicActionD(action, PParams, ml, mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<WilsonImplF>* paction_f_p1 = createPeriodicActionF(action, PParams, ml, mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  CayleyFermion5D<WilsonImplD>* paction_d_p2(paction_d_p1);
  CayleyFermion5D<WilsonImplF>* paction_f_p2(paction_f_p1);
#endif


  
  LatticeSCFmatrixD src1_GP = fourMomentumVolumeSource<GparityWilsonImplD>(minusMom(p1_GP), UGridD); //we use exp(+i) Fourier transform conventions for source \bar\psi whereas the code uses exp(-i)
  LatticeSCFmatrixD src2_GP = fourMomentumVolumeSource<GparityWilsonImplD>(minusMom(p2_GP), UGridD);
  
  MixedCGargs cg_args;
  LatticeSCFmatrixD prop1_GP = mixedPrecInvert(src1_GP, *action_d, *action_f, cg_args);
  LatticeSCFmatrixD prop2_GP = mixedPrecInvert(src2_GP, *action_d, *action_f, cg_args);
  
  SCFmatrixD prop1_GP_sinkFT = sinkTransform(p1_GP, prop1_GP);
  SCFmatrixD prop2_GP_sinkFT = sinkTransform(p2_GP, prop2_GP);

  SCFmatrixD prop1_GP_sinkFT_inv = invert( prop1_GP_sinkFT );
  SCFmatrixD prop2_GP_sinkFT_inv = invert( prop2_GP_sinkFT );

  LatticeSCFmatrixD D_p1_in = mulF11leftright( LatticeSCFmatrixD(prop1_GP * prop1_GP_sinkFT_inv ) );
  LatticeSCFmatrixD D_p2_in = mulF11leftright( LatticeSCFmatrixD(prop2_GP * prop2_GP_sinkFT_inv ) );
  
  Gamma gamma5(Gamma::Algebra::Gamma5);

  LatticeSCFmatrixD D_p1_out = gamma5 * adj(D_p1_in) * gamma5;

  std::cout << "|D_p1_out|^2=" << norm2(D_p1_out) << " |D_p2_in|^2=" << norm2(D_p2_in) << std::endl;

  //Do scalar vertex
  LatticeComplexD vphase = phaseFieldFour(minusMom(p1_GP), UGridD) * phaseFieldFour(p2_GP, UGridD); //exp(i(p1-p2)x)

  LatticeSCFmatrixD Ax = vphase * D_p1_out * (GammaOp * D_p2_in);
  SCFmatrixD As = GammaOp * sum(Ax);
  ComplexD A = trace(As);

  //Periodic
  LatticeSpinColourMatrixD src1_P = fourMomentumVolumeSource<WilsonImplD>(minusMom(p1_P), UGridD);
  LatticeSpinColourMatrixD src2_P = fourMomentumVolumeSource<WilsonImplD>(minusMom(p2_P), UGridD);
  
  LatticeSpinColourMatrixD prop1_P = mixedPrecInvert(src1_P, *paction_d_p1, *paction_f_p1, cg_args);
  LatticeSpinColourMatrixD prop2_P = mixedPrecInvert(src2_P, *paction_d_p2, *paction_f_p2, cg_args);
  
  SpinColourMatrixD prop1_P_sinkFT = sinkTransform(p1_P, prop1_P);
  SpinColourMatrixD prop2_P_sinkFT = sinkTransform(p2_P, prop2_P);

  SpinColourMatrixD prop1_P_sinkFT_inv = invert( prop1_P_sinkFT );
  SpinColourMatrixD prop2_P_sinkFT_inv = invert( prop2_P_sinkFT );

  LatticeSpinColourMatrixD Gamp_p1_in_P = prop1_P * prop1_P_sinkFT_inv;
  LatticeSpinColourMatrixD Gamp_p2_in_P = prop2_P * prop2_P_sinkFT_inv;
  
  LatticeSpinColourMatrixD Gamp_p1_out_P = gamma5 * adj(Gamp_p1_in_P) * gamma5;
  
  //Do scalar vertex
  LatticeComplexD vphase_P = phaseFieldFour(minusMom(p1_P), UGridD) * phaseFieldFour(p2_P, UGridD);

  LatticeSpinColourMatrixD Ax_P = vphase_P * Gamp_p1_out_P * (GammaOp * Gamp_p2_in_P);
  SpinColourMatrixD As_P = GammaOp * sum(Ax_P);
  ComplexD A_P = trace( As_P);

  std::cout << "Gparity: " << A << " Periodic: " << A_P << std::endl;

  return 0;
}
