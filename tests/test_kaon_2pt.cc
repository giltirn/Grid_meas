#include<common.h>

using namespace Grid;
using namespace GridMeas;
 

int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  int Ls = 12;
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];
  size_t V4d = latt[0]*latt[1]*latt[2]*latt[3];
  size_t V3d = latt[0]*latt[1]*latt[2];

  auto UGridD   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  auto UrbGridD = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridD);
  auto FGridD     = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridD);
  auto FrbGridD   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridD);

  auto UGridF   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  auto UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls, UGridF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGridF);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(UGridD); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD U_d(UGridD);
  LatticeGaugeFieldF U_f(UGridF);

  //Use a unit gauge config to test various properties
  SU<Nc>::ColdConfiguration(U_d);
  precisionChange(U_f, U_d);
  
  //Setup BCs
  std::vector<int> GparityDirs = {1,0,0};

  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction

  CayleyFermion5D<GparityWilsonImplD>* action_d = createActionD(ActionType::DWF, Params, 0.01, 1., U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<GparityWilsonImplF>* action_f = createActionF(ActionType::DWF, Params, 0.01, 1., U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  CayleyFermion5D<GparityWilsonImplD>* action_s_d = createActionD(ActionType::DWF, Params, 0.032, 1., U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<GparityWilsonImplF>* action_s_f = createActionF(ActionType::DWF, Params, 0.032, 1., U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  std::vector<int> p(3,0);
  for(int i=0;i<3;i++)
    p[i] = GparityDirs[i] ? 1 : 0;
  std::vector<double> p_phys = getPhysicalMomentum(p);

  std::cout << "p=" << p << "  p_phys=" << p_phys << std::endl;
  
  LatticeComplexD p_src_phase_field = phaseField(p_phys, UGridD); //exp(-i \vec p \cdot \vec x)
  p_src_phase_field = conjugate(p_src_phase_field); //exp(+i \vec p \cdot \vec x)

  //Test translational invariance of 2pt
  LatticeSCFmatrixD src_t1 = wallSource(1, UGridD);
  src_t1 = src_t1 * p_src_phase_field;

  LatticeSCFmatrixD src_t4 = wallSource(4, UGridD);
  src_t4 = src_t4 * p_src_phase_field;
  
  LatticeSCFmatrixD R_t1(UGridD);
  R_t1 = mixedPrecInvert(src_t1, *action_d, *action_f, 1e-10, 1e-6);

  LatticeSCFmatrixD R_h_t1(UGridD);
  R_h_t1 = mixedPrecInvert(src_t1, *action_s_d, *action_s_f, 1e-10, 1e-6);

  {
    LatticeSCFmatrixD test = R_t1 * adj(R_h_t1);
    std::vector<SCFmatrixD> C;
    sliceSum(test, C, 3);
    std::cout << GridLogMessage << "X(t) = \\sum_\\vec x tr(R(x)R_h^dag(x)) for prop with t0=" << 1 << std::endl;
    for(int t=0;t<Lt;t++){
      std::cout << t << " " << trace(C[t]) << std::endl;
    }
  }

  LatticeSCFmatrixD R_t4(UGridD);
  R_t4 = mixedPrecInvert(src_t4, *action_d, *action_f, 1e-10, 1e-6);

  LatticeSCFmatrixD R_h_t4(UGridD);
  R_h_t4 = mixedPrecInvert(src_t4, *action_s_d, *action_s_f, 1e-10, 1e-6);


  {
    LatticeSCFmatrixD test = R_t4 * adj(R_h_t4);
    std::vector<SCFmatrixD> C;
    sliceSum(test, C, 3);
    std::cout << GridLogMessage << "X(t) = \\sum_\\vec x tr(R(x)R_h^dag(x)) for prop with t0=" << 4 << std::endl;
    for(int t=0;t<Lt;t++){
      std::cout << t << " " << trace(C[t]) << std::endl;
    }
  }

  std::vector<RealD> Cts_t1 = momWallSourceKaonCorrelator(p, 1, R_t1, R_h_t1);
  std::vector<RealD> Cts_t4 = momWallSourceKaonCorrelator(p, 4, R_t4, R_h_t4);

  assert(Cts_t1.size() == Lt);
  assert(Cts_t4.size() == Lt);
  
  for(int i=0;i<Lt;i++){
    std::cout << "C(" << i << ", t0=" << 1 << ")=" << Cts_t1[i] << "    C(" << i << ", t0=" << 4 << ")=" << Cts_t4[i] << std::endl;    
    assert(fabs( Cts_t1[i] ) > 1e-2 ); //cannot just be zero within noise (this will happen if you mess up the momentum projection!)
    assert(fabs( Cts_t1[i] - Cts_t4[i] ) < 1e-6);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
