#include<common/propagator_invert.h>
#include<common/utils.h>
#include<common/grids.h>
#include<common/action.h>
#include<common/sources.h>

using namespace Grid;
using namespace GridMeas;

 

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  int Ls = 12;
  Coordinate latt = GridDefaultLatt();
  Coordinate proc_full = GridDefaultMpi();

  //We need at least 2 ranks for this to be a meaningful test
  int ranks = 1; for(int i=0;i<Nd;i++) ranks *= proc_full[i];
  assert(ranks > 1);
  
  //Basic setup
  Grids fullGridsD = makeDoublePrecGrids(Ls, latt);
  Grids fullGridsF = makeSinglePrecGrids(Ls, latt);

  std::vector<int> seeds4({0,1,2,3});
  GridParallelRNG pRNG(fullGridsD.UGrid); //4D!
  pRNG.SeedFixedIntegers(seeds4);
  
  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  double mobius_scale = 2.;
  double ml = 0.01;
  ActionType action = ActionType::Mobius;
  
  assert(Nd == 4);
  std::vector<int> dirs4(Nd);
  std::vector<int> GparityDirs = {1,1,1};
  for(int i=0;i<3;i++) dirs4[i] = GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction
  std::cout << GridLogMessage << "Fermion BCs: " << Params.twists << std::endl;

  LatticeGaugeFieldD U_d(fullGridsD.UGrid);
  SU<Nc>::HotConfiguration(pRNG, U_d);

  LatticeGaugeFieldF U_f(fullGridsF.UGrid);
  precisionChange(U_f, U_d);
 
  CayleyFermion5D<GparityWilsonImplD>* action_d = createActionD(action, Params, ml, mobius_scale, U_d, fullGridsD);
  CayleyFermion5D<GparityWilsonImplF>* action_f = createActionF(action, Params, ml, mobius_scale, U_f, fullGridsF);

  //For convenient have 1-rank subgrids
  std::cout << GridLogMessage << "Generating split grids" << std::endl;
  Coordinate proc_sub(std::vector<int>({1,1,1,1}));
  Grids subGridsD = makeSplitGrids(fullGridsD, proc_sub);
  Grids subGridsF = makeSplitGrids(fullGridsF, proc_sub);

  LatticeGaugeFieldD U_sub_d(subGridsD.UGrid);
  Grid_split(U_d,U_sub_d);
  LatticeGaugeFieldF U_sub_f(subGridsF.UGrid);
  precisionChange(U_sub_f,U_sub_d);
  
  CayleyFermion5D<GparityWilsonImplD>* action_sub_d = createActionD(action, Params, ml, mobius_scale, U_sub_d, subGridsD);
  CayleyFermion5D<GparityWilsonImplF>* action_sub_f = createActionF(action, Params, ml, mobius_scale, U_sub_f, subGridsF);
  
  //Create 2 sources
  std::vector<LatticeSCFmatrixD> src(2, fullGridsD.UGrid);
  src[0]= Z2wallSource(0, pRNG, fullGridsD.UGrid);
  src[1]= Z2wallSource(0, pRNG, fullGridsD.UGrid);
    
  std::vector<LatticeSCFmatrixD> sol_basic(2, fullGridsD.UGrid);
  std::vector<LatticeSCFmatrixD> sol_split(2, fullGridsD.UGrid);
  
  //Do regular solves
  std::cout << GridLogMessage << "Performing regular inversion" << std::endl;
  for(int i=0;i<2;i++)
    sol_basic[i] = mixedPrecInvert(src[i], *action_d, *action_f, 1e-8, 1e-5);

  //Do split solves
  std::cout << GridLogMessage << "Performing split grid inversion" << std::endl;
  splitGridMixedPrecInvert(sol_split, src, *action_d, *action_sub_d, *action_sub_f, 1e-8, 1e-5);

  //Compare
  for(int i=0;i<2;i++){
    LatticeSCFmatrixD diff = sol_basic[i] - sol_split[i];
    std::cout << GridLogMessage << norm2(sol_basic[i]) << " " << norm2(sol_split[i]) << " diff " << norm2(diff) << std::endl;
  }
  
  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}

	  
