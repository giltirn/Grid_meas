#include<common/propagator_invert.h>
#include<common/propagator_invert_Xconj.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/momentum.h>
#include<common/sources.h>
#include<common/pion_2pt.h>
#include<common/ps_singlet_2pt.h>

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

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD U_d(GridsD.UGrid);
  LatticeGaugeFieldF U_f(GridsF.UGrid);
  SU<Nc>::HotConfiguration(pRNG,U_d);
  precisionChange(U_f,U_d);
  
  std::vector<Integer> GparityDirs = {1,1,1};
  GparityWilsonImplD::ImplParams Params = setupActionParams(GparityDirs);
  
  Actions actions(ActionType::Mobius, Params, 0.01, 2.0, U_d, GridsD, U_f, GridsF);

  for(int pm=0;pm<2;pm++){   
    //Choose a G-parity momentum
    std::vector<int> p1(3,0);
    for(int i=0;i<3;i++){
      p1[i] = GparityDirs[i] ? 1 : 0;
      if(pm) p1[i] = -p1[i];
    }

    std::cout << "Quark momentum " << momstr_human(p1) << std::endl;

    std::vector<double> p1_phys = getPhysicalMomentum(p1);

    //Compute a propagator with momentum sources
    int t0 = 0;
    MixedCGargs cg_args;
    LatticeSCFmatrixD src_wall = momentumWallSource(p1_phys, t0, GridsD.UGrid);
    LatticeSCFmatrixD prop_wall = mixedPrecInvert(src_wall, *actions.action_d, *actions.action_f, cg_args);

    LatticeSCFmatrixD src_gpcos = gparityCosineWallSource(p1, t0, GridsD.UGrid);
    LatticeSCFmatrixD prop_gpcos = mixedPrecInvertXconj(src_gpcos, *actions.xconj_action_d, *actions.xconj_action_f, cg_args);

    std::vector<RealD> pion_wall =  momWallSourcePionCorrelator(p1,p1,t0,prop_wall,prop_wall);
    std::vector<RealD> pion_gpcos =  momWallSourcePionCorrelator(p1,p1,t0,prop_gpcos,prop_gpcos);

    std::vector<RealD> ps_singlet_wall =  momWallSourcePSsingletCorrelator(p1, t0, prop_wall);
    std::vector<RealD> ps_singlet_gpcos =  momWallSourcePSsingletCorrelator(p1, t0, prop_gpcos);

    std::cout << "Test pion 2pt, quark momentum " << momstr_human(p1) << std::endl;
    for(int t=0;t<latt[3];t++){
      RealD reldiff = 2*(pion_wall[t] - pion_gpcos[t])/(pion_wall[t] + pion_gpcos[t]);
      std::cout << t << " wall: " << pion_wall[t] << " gpcos: " << pion_gpcos[t] << " reldiff: " << reldiff << std::endl;
      assert(fabs(reldiff) < 1e-6 );
    }

    std::cout << "Test pseudoscalar single 2pt, quark momentum " << momstr_human(p1) << std::endl;
    for(int t=0;t<latt[3];t++){
      RealD reldiff = 2*(ps_singlet_wall[t] - ps_singlet_gpcos[t])/(ps_singlet_wall[t] + ps_singlet_gpcos[t]);
      std::cout << t << " wall: " << ps_singlet_wall[t] << " gpcos: " << ps_singlet_gpcos[t] << " reldiff: " << reldiff << std::endl;
      assert(fabs(reldiff) < 1e-6 );
    }
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}

	  
