#include<common/propagator_invert_field_Xconj.h>
#include<common/propagator_invert_field.h>
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
  
  FermionFieldD src(GridsD.UGrid);
  gaussian(pRNG,src); //no restrictions on src

  FermionFieldD sol_GP(GridsD.UGrid), sol_Xconj(GridsD.UGrid);

  mixedPrecInvertField(sol_GP, src, *actions.action_d, *actions.action_f, 1e-8, 1e-5, (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
  mixedPrecInvertFieldXconj(sol_Xconj, src, *actions.xconj_action_d, *actions.xconj_action_f, 1e-8, 1e-5, (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);

  FermionFieldD diff = sol_Xconj - sol_GP;
  std::cout << "Norm of difference (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
