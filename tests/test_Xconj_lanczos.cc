#include<common/propagator_invert.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/lanczos.h>

using namespace Grid;
using namespace GridMeas;

  
  

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  int Ls = 8;
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];
  size_t V4d = latt[0]*latt[1]*latt[2]*latt[3];
  size_t V3d = latt[0]*latt[1]*latt[2];

  Grids GridsD = makeDoublePrecGrids(Ls, latt);
  Grids GridsF = makeSinglePrecGrids(Ls, latt);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(GridsD.UGrid); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridParallelRNG pRNG5(GridsD.FGrid);  
  pRNG5.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD Ud(GridsD.UGrid);
  LatticeGaugeFieldF Uf(GridsF.UGrid);
  SU<Nc>::HotConfiguration(pRNG, Ud);
  precisionChange(Uf,Ud);
  
  std::vector<int> dirs4({1,1,1,0});
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  Params.twists = Coordinate(std::vector<int>({1,1,1,1})); //APBC in t direction

  Actions actions(ActionType::Mobius, Params, 0.01, 2.0, Ud, GridsD, Uf, GridsF);
  
  LanczosParameters params;
  
  int Nstop = 30;
  int Nk = 30;
  int Np = 5;
  int ord = 80;
  RealD lo = 1.5;
  RealD hi = 88.0;
  params.alpha = sqrt(hi);
  params.beta = sqrt(lo);
  params.n_stop = Nstop;
  params.n_want = Nk;
  params.n_use = Nk+Np;
  params.ord = ord;
  params.tolerance = 1e-8;

  std::vector<RealD> eval_GP;
  std::vector<FermionFieldD> evec_GP;

  computeEigenvalues(eval_GP, evec_GP, params, GridsD.FGrid, GridsD.FrbGrid, Ud, *actions.action_d, pRNG5);


  std::vector<RealD> eval_Xconj;
  std::vector<FermionField1fD> evec_Xconj;

  computeEigenvalues(eval_Xconj, evec_Xconj, params, GridsD.FGrid, GridsD.FrbGrid, Ud, *actions.xconj_action_d, pRNG5);

  assert(eval_GP.size() == eval_Xconj.size());

  std::cout << "Comparing evals: " << std::endl;
  for(int i=0;i<eval_GP.size();i++){
    std::cout << eval_Xconj[i] << " " << eval_GP[i] << " diff: " << eval_Xconj[i] - eval_GP[i] << std::endl;
  }

  Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  Gamma X = C*g5;

  std::cout << "Comparing eigenvectors: " << std::endl;
  for(int i=0;i<eval_GP.size();i++){
    //We need to phase-rotate the regular evecs into X-conjugate vectors
    FermionField1fD v0 = PeekIndex<GparityFlavourIndex>(evec_GP[i],0);
    FermionField1fD v1 = PeekIndex<GparityFlavourIndex>(evec_GP[i],1);
    FermionField1fD Xv1star = X*conjugate(v1);
    ComplexD z = innerProduct(v0, Xv1star);
    ComplexD alpha = ComplexD(0.5)/z;

    FermionFieldD w = 1./sqrt(alpha) * evec_GP[i];
    //w should be X-conjugate; check
    FermionField1fD w0 = PeekIndex<GparityFlavourIndex>(w,0);
    FermionField1fD w1 = PeekIndex<GparityFlavourIndex>(w,1);
    FermionField1fD tmp = w1 + X*conjugate(w0);
    std::cout << "Converted regular evec to X-conjugate: check (expect 0): " << norm2(tmp) << std::endl;

    tmp = w0 - evec_Xconj[i];
    FermionField1fD tmp2 = w0 + evec_Xconj[i];
    std::cout << "Evec " << i << " difference: " << norm2(tmp) << " sum: " << norm2(tmp2) << std::endl;
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
