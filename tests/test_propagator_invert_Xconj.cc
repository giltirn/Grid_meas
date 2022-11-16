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
  
  //Should work with any spin-color source phi providing its flavor structure is diag(phi, phi*)  and [phi,X] = 0
  static Gamma X = Xmatrix();

  LatticeSpinColourMatrixD phi(GridsD.UGrid);
  gaussian(pRNG, phi);
  phi = X*phi + phi*X;

  LatticeSCFmatrixD src(GridsD.UGrid);
  src = Zero();
  PokeIndex<GparityFlavourIndex>(src, phi, 0,0);
  PokeIndex<GparityFlavourIndex>(src, conjugate(phi), 1,1);
  
  //Test left and right multiplication by H
  {
    LatticeSCFmatrixD r1 = mulHdagLeft(mulHLeft(src));
    LatticeSCFmatrixD r2 = mulHLeft(mulHdagLeft(src));
    LatticeSCFmatrixD r3 = mulHdagRight(mulHRight(src));
    LatticeSCFmatrixD r4 = mulHRight(mulHdagRight(src));
    LatticeSCFmatrixD tmp(src.Grid());
    tmp = r1 - src;
    RealD diff1 = norm2(tmp);
    tmp = r2 - src;
    RealD diff2 = norm2(tmp);
    tmp = r3 - src;
    RealD diff3 = norm2(tmp);
    tmp = r4 - src;
    RealD diff4 = norm2(tmp);

    std::cout << "Test H and H^-1 left and right multiplication (expect 0): " << diff1 << " " << diff2 << " " << diff3 << " " << diff4 << std::endl;
    assert(diff1 < 1e-8);
    assert(diff2 < 1e-8);
    assert(diff3 < 1e-8);
    assert(diff4 < 1e-8);
  }
  MixedCGargs cg_args;

  LatticeSCFmatrixD sol_GP(GridsD.UGrid), midsol_GP(GridsD.UGrid);
  mixedPrecInvertGen(sol_GP, midsol_GP, src, *actions.action_d, *actions.action_f, cg_args, true, 
		     (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);

  LatticeSCFmatrixD sol_Xconj(GridsD.UGrid), midsol_Xconj(GridsD.UGrid);
  mixedPrecInvertGenXconj(sol_Xconj, midsol_Xconj, src, *actions.xconj_action_d, *actions.xconj_action_f, cg_args, true, 
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
