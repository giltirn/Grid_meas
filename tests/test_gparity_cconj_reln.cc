#include<common/propagator_invert.h>
#include<common/sources.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/momentum.h>

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

  assert(Nd == 4);
  std::vector<int> dirs4({1,0,0,0});
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  Params.twists = Coordinate(std::vector<int>({1,0,0,1})); //APBC in t direction

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(UGridD); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD U_d(UGridD);
  LatticeGaugeFieldF U_f(UGridF);
  gaussian(pRNG, U_d);
  precisionChange(U_f, U_d);
  
  ActionType action = ActionType::Mobius;
  RealD ml = 0.01;
  RealD mobius_scale = 2.0;

  CayleyFermion5D<GparityWilsonImplD>* action_d = createActionD(action, Params, ml, mobius_scale, U_d, *FGridD, *FrbGridD, *UGridD, *UrbGridD);
  CayleyFermion5D<GparityWilsonImplF>* action_f = createActionF(action, Params, ml, mobius_scale, U_f, *FGridF, *FrbGridF, *UGridF, *UrbGridF);

  GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
  GparityFlavour sigma2(GparityFlavour::Algebra::SigmaY);
  Gamma gamma5(Gamma::Algebra::Gamma5);
  Gamma gamma2(Gamma::Algebra::GammaY);
  Gamma mgamma2(Gamma::Algebra::MinusGammaY);
  Gamma gamma4(Gamma::Algebra::MinusGammaT);
  

  {
    LatticeSCFmatrixD rsrc = Z2wallSource(1, pRNG, UGridD); //source needs to be real and have unit spin and flavor structure
    LatticeSCFmatrixD prop = mixedPrecInvert(rsrc, *action_d, *action_f, 1e-8, 1e-5);
    
    //G* = g5 Cinv sigma2 G sigma2 C g5
    //C = -g2 g4

    
    //C=-\gamma^2 \gamma^4
    //C^{-1} = -C
    LatticeSCFmatrixD rhs = gamma5 * ( gamma2 * ( gamma4 * (sigma2 * prop * sigma2) * mgamma2 ) * gamma4 ) * gamma5;
    LatticeSCFmatrixD lhs = conjugate(prop);
    LatticeSCFmatrixD diff = rhs - lhs;
    
    RealD norm2_diff = norm2(diff);
    std::cout << "CConj relation on real Z2wall propagator: lhs=" << norm2(lhs) << " rhs=" << norm2(rhs) << " diff=" << norm2_diff << " expect 0" << std::endl;
    assert(norm2_diff < 1e-10);
  }

  {
    //Test the following:
    //R(\vec x, t; \vec p_1) = \sum_{\vec y_1} e^{i \vec p_1 \cdot \vec y_1} G(\vec x, t; \vec y_1, t_0)
    //S(\vec x, t; \vec p_2) = \sum_{\vec y_2} e^{i \vec p_2 \cdot \vec y_2} \gamma^5 G^\dagger(\vec x, t; \vec y_2, t_0) \gamma^5
    //             = \sum_{\vec y_2} e^{i \vec p_2 \cdot \vec y_2} \gamma^5 [ \gamma^5 C^{-1} \sigma_2 G^T(\vec x, t; \vec y_2, t_0) C\gamma^5 \sigma_2] \gamma^5
    //             = \sum_{\vec y_2} e^{i \vec p_2 \cdot \vec y_2} C^{-1} \sigma_2 G^T(\vec x, t; \vec y_2, t_0) C \sigma_2
    //             = C^{-1} \sigma_2 [ \sum_{\vec y_2} e^{i \vec p_2 \cdot \vec y_2} G(\vec x, t; \vec y_2, t_0) ]^T C \sigma_2
    //             = C^{-1} \sigma_2 R(\vec x, t; \vec p_2)^T C \sigma_2
    
    LatticeSCFmatrixD z2_src = Z2wallSource(1, pRNG, UGridD); //source needs to be real and have unit spin and flavor structure
    std::vector<int> p = {1,0,0};
    std::vector<double> pphys = getPhysicalMomentum(p);

    std::vector<int> mp = {-1,0,0};
    std::vector<double> mpphys = getPhysicalMomentum(mp);

    //exp(i \vec p \cdot \vec x)
    LatticeComplexD p_phase_field = phaseField(pphys, UGridD); 
    conjugate(p_phase_field);

    //exp(-i \vec p \cdot \vec x)
    LatticeComplexD mp_phase_field = phaseField(mpphys, UGridD); 
    conjugate(mp_phase_field);
    
    //Compute R
    LatticeSCFmatrixD Rpsrc = z2_src * p_phase_field;
    
    LatticeSCFmatrixD Rp = mixedPrecInvert(Rpsrc, *action_d, *action_f, 1e-8, 1e-5);
    
    //Construct S the hard way
    //S(\vec x, t; \vec p_2) = \gamma^5 [ \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2}  G(\vec x, t; \vec y_2, t_0) ] ^\dagger\gamma^5

    LatticeSCFmatrixD Spsrc = z2_src * mp_phase_field;
    LatticeSCFmatrixD Sp = mixedPrecInvert(Spsrc, *action_d, *action_f, 1e-8, 1e-5);
    Sp = adj(Sp); //-p -> p
    Sp = gamma5 * Sp * gamma5;

    //Construct S from R
    LatticeSCFmatrixD Srecon = gamma2 * ( gamma4 * ( sigma2 * transpose(Rp) * mgamma2 ) * gamma4 ) * sigma2;

    LatticeSCFmatrixD diff = Sp - Srecon;
    RealD n = norm2(diff);

    std::cout << "CConj relation relate S and R in momentum wall sources : lhs=" << norm2(Sp) << " rhs=" << norm2(Srecon) << " diff=" << n << " expect 0" << std::endl;
    assert(n < 1e-7);
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}

