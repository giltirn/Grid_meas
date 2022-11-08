#include<common/propagator_invert.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/momentum.h>
#include<common/sources.h>
#include<common/pion_2pt.h>

using namespace Grid;
using namespace GridMeas;

//Use wrong projector for testing
std::vector<RealD> momWallSourcePionCorrelatorOppProj(const std::vector<int> &p1, const std::vector<int> &p2, int t0, const LatticeSCFmatrixD& R_p1, const LatticeSCFmatrixD& R_p2){
std::cout << GridLogMessage << "Starting momentum wall source pion correlator with t0=" << t0 << std::endl;
 GridBase* UGrid = R_p1.Grid();
 std::vector<int> p = { p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2] };
 std::vector<double> pphys = getPhysicalMomentum(p);
 std::cout << GridLogMessage << "Total momentum in base units = " << p << " and physical units = " << pphys << std::endl;
 LatticeComplexD snk_phase_field = conjugate(phaseField(pphys, R_p1.Grid())); //exp(i \vec p \cdot \vec x)
  
 GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
 GparityFlavour sigma2(GparityFlavour::Algebra::SigmaY);
 Gamma gamma5(Gamma::Algebra::Gamma5);
 Gamma gamma2(Gamma::Algebra::GammaY);
 Gamma mgamma2(Gamma::Algebra::MinusGammaY);
 Gamma gamma4(Gamma::Algebra::MinusGammaT);
  
 //C=-\gamma^2 \gamma^4
 //C^{-1} = -C
 LatticeSCFmatrixD S_p2 = gamma2 * ( gamma4 * ( sigma2 * transpose(R_p2) * mgamma2 ) * gamma4 )* sigma2;

 GparityFlavour::Algebra proj = getProjector(p2);
 proj = ( proj == GparityFlavour::Algebra::ProjPlus ? GparityFlavour::Algebra::ProjMinus : GparityFlavour::Algebra::ProjPlus );
 GparityFlavour Xi(proj);

 std::cout << GridLogMessage << "Projector for p2=" << p2 << " : " << proj << std::endl;

 LatticeSCFmatrixD sqb = snk_phase_field * ( sigma3* ( gamma5 * R_p1 * gamma5 ) * sigma3 ) * Xi * S_p2;
  
 std::vector<SCFmatrixD> Ctm;
 sliceSum(sqb, Ctm, 3);
  
 int Lt = Ctm.size();
 std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

 for(int t=0;t<Lt;t++)
   std::cout << GridLogMessage << "C(t)=" << 0.5*trace(Ctm[t]) << std::endl;
    
 std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
 for(int t=0;t<Lt;t++)
   Ctx[t] = 0.5 * real( trace(Ctm[t]) );

 std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
 for(int tx=0;tx<Lt;tx++){
   int t = ( tx - t0 + Lt ) % Lt;
   Ct[t] = Ctx[tx];
 }
 return Ct;
}




int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  bool coswall = false;
  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "--coswall") coswall = true;
  }

  int Ls = 12;
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];
  size_t V4d = latt[0]*latt[1]*latt[2]*latt[3];
  size_t V3d = latt[0]*latt[1]*latt[2];

  Grids GridsD = makeDoublePrecGrids(Ls, latt);
  Grids GridsF = makeSinglePrecGrids(Ls, latt);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(GridsD.UGrid); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  //We are testing boundary translational covariance so a unit gauge is sufficient
  LatticeGaugeFieldD U_d(GridsD.UGrid);
  LatticeGaugeFieldF U_f(GridsF.UGrid);
  SU<Nc>::ColdConfiguration(U_d);
  precisionChange(U_f,U_d);
  
  std::vector<Integer> GparityDirs = {1,1,1};
  GparityWilsonImplD::ImplParams Params = setupActionParams(GparityDirs);
  
  Actions actions(ActionType::Mobius, Params, 0.01, 2.0, U_d, GridsD, U_f, GridsF);

  //Choose a G-parity momentum
  std::vector<int> p1(3,0);
  for(int i=0;i<3;i++)
    p1[i] = GparityDirs[i] ? 1 : 0;
  std::vector<double> p1_phys = getPhysicalMomentum(p1);

  std::vector<int> mp1({-p1[0],-p1[1],-p1[2]});

  //Compute a propagator with momentum sources
  int t0 = 0;
  
  LatticeSCFmatrixD eta = coswall ? cosineWallSource(p1_phys, t0, GridsD.UGrid) : momentumWallSource(p1_phys, t0, GridsD.UGrid);
  LatticeSCFmatrixD prop = mixedPrecInvert(eta, *actions.action_d, *actions.action_f, 1e-8, 1e-5);

  std::vector<RealD> pion_rightproj_rightmom =  momWallSourcePionCorrelator(p1,p1,t0,prop,prop);
  std::vector<RealD> pion_wrongproj_rightmom =  momWallSourcePionCorrelatorOppProj(p1,p1,t0,prop,prop);
  std::vector<RealD> pion_wrongproj_wrongmom =  momWallSourcePionCorrelator(mp1,mp1,t0,prop,prop);
  std::vector<RealD> pion_rightproj_wrongmom =  momWallSourcePionCorrelatorOppProj(mp1,mp1,t0,prop,prop);

  if(coswall) std::cout << "For cosine sources both +/-p sink momenta are valid, but if we use the wrong projector we should get a smaller projection" << std::endl;
  else std::cout << "The right projector will project only onto the right momentum +p, whereas the wrong projector will project onto +/-p" << std::endl;

  for(int t=0;t<latt[3];t++){
    std::cout << t << " rightproj,rightmom: " << pion_rightproj_rightmom[t] << " wrongproj,rightmom: " << pion_wrongproj_rightmom[t]
	      << " rightproj,wrongmom: " << pion_rightproj_wrongmom[t] << " wrongproj,wrongmom: " << pion_wrongproj_wrongmom[t] << std::endl; 
    if(!coswall){
      assert(fabs(pion_rightproj_wrongmom[t]) < 1e-6);
      assert(fabs(pion_wrongproj_wrongmom[t]) > 1e-4);
    }
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}

	  
