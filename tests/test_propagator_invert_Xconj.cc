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
  
  //Should work with any matrix of the form
  //|  A       B   |
  //| XB*X   -XA*X |

  LatticePropagator A(GridsD.UGrid);
  gaussian(pRNG, A);

  LatticePropagator B(GridsD.UGrid);
  gaussian(pRNG, B);
  
  LatticePropagator C = Xmatrix()*(conjugate(B)*Xmatrix());
  LatticePropagator D = -(Xmatrix()*(conjugate(A)*Xmatrix()));
  
  LatticeSCFmatrixD src(GridsD.UGrid);
  src = Zero();

  PokeIndex<GparityFlavourIndex>(src, A, 0,0);
  PokeIndex<GparityFlavourIndex>(src, B, 0,1);
  PokeIndex<GparityFlavourIndex>(src, C, 1,0);
  PokeIndex<GparityFlavourIndex>(src, D, 1,1);
 
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
  //Test src*H has X-conjugate columns
  {
    std::cout << "Test src*H has X-conjugate columns:" << std::endl;
    LatticeSCFmatrixD r = mulHRight(src);
    typedef columnOps<FermionFieldD> cop;
    for(int c=0;c<24;c++){
      FermionFieldD col = cop::extractColumn(r, c);
      FermionField1fD col_0 = PeekIndex<GparityFlavourIndex>(col,0);
      FermionField1fD col_1 = PeekIndex<GparityFlavourIndex>(col,1);
      FermionField1fD tmp = col_1 + Xmatrix()*conjugate(col_0);
      std::cout << c << " (expect 0): " << norm2(tmp) << std::endl;
      assert(norm2(tmp) < 1e-10);
    }
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
