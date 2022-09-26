#include<common/propagator_invert_Xconj.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/sources.h>

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
  
  //LatticeSCFmatrixD src = wallSource(0, GridsD.UGrid);

  //Should work with any spin-color source phi providing its flavor structure is diag(phi, phi*)  and [phi,X] = 0
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;

  LatticeSpinColourMatrixD phi(GridsD.UGrid);
  gaussian(pRNG, phi);
  phi = X*phi + phi*X;

  LatticeSCFmatrixD src(GridsD.UGrid);
  src = Zero();
  PokeIndex<GparityFlavourIndex>(src, phi, 0,0);
  PokeIndex<GparityFlavourIndex>(src, conjugate(phi), 1,1);
  
  LatticeSCFmatrixD sol_GP(GridsD.UGrid), midsol_GP(GridsD.UGrid);
  mixedPrecInvertGen(sol_GP, midsol_GP, src, *actions.action_d, *actions.action_f, 1e-8, 1e-5, true, 
		     (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);

  LatticeSCFmatrixD sol_Xconj(GridsD.UGrid), midsol_Xconj(GridsD.UGrid);
  mixedPrecInvertGenXconj(sol_Xconj, midsol_Xconj, src, *actions.xconj_action_d, *actions.xconj_action_f, 1e-8, 1e-5, true, 
  			  (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);

  LatticeSCFmatrixD diff = sol_Xconj - sol_GP;
  std::cout << "Prop: Norm of difference (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);
  diff = midsol_Xconj - midsol_GP;
  std::cout << "Midprop: Norm of difference (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);

  



  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
