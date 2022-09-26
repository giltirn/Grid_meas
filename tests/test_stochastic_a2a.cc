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

  Coordinate site({1,1,1,1}), other_site({0,1,0,1});
  FermionFieldD point_src(GridsD.UGrid);
  point_src = Zero();
  
  typename SiteSpinorD::scalar_object pt;
  pt(0)(0)(0) = 1.0;
  pokeSite(pt,point_src,site);

  FermionFieldD point_sol(GridsD.UGrid);
  mixedPrecInvertField(point_sol, point_src, *actions.action_d, *actions.action_f, 1e-8, 1e-5, (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);

  FermionFieldD rsum(GridsD.UGrid), wsum(GridsD.UGrid);
  rsum = Zero();
  wsum = Zero();

  int hits = 20000;
  for(int h=0;h<hits;h++){
    FermionFieldD rnd_src = randomGaussianVolumeSource(pRNG, GridsD.UGrid);
    FermionFieldD rnd_sol(GridsD.UGrid);
    mixedPrecInvertField(rnd_sol, rnd_src, *actions.action_d, *actions.action_f, 1e-8, 1e-5, (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
    //rnd_sol = rnd_src;

    typename SiteSpinorD::scalar_object src_p = conjugate(peekSite(rnd_src, site));
    
    rsum = rsum + rnd_sol * src_p(0)(0)(0); //v(x) w*(y)
    wsum = wsum + rnd_src * src_p(0)(0)(0); //w(x) w*(y)  ~ delta_x,y -> check gives 1 for source site
    

    FermionFieldD ravg = rsum * RealD(1./(h+1));
    FermionFieldD diff = point_sol - ravg;
    typename SiteSpinorD::scalar_object  v = peekSite(ravg, other_site);
    typename SiteSpinorD::scalar_object  ve = peekSite(point_sol, other_site);

    FermionFieldD wavg = wsum * RealD(1./(h+1));
    typename SiteSpinorD::scalar_object  w = peekSite(wavg, site);
    typename SiteSpinorD::scalar_object  we = peekSite(wavg, other_site);

    std::cout << "Hit " << h << " running norm diff " << norm2(diff) << " a sink site val " << v(0)(0)(0) << " expect " << ve(0)(0)(0) << " ratio " <<  v(0)(0)(0)/ve(0)(0)(0) 
	      << ". Check src delta norm, src site : " <<  w(0)(0)(0) << " (expect 1), some other site : " << we(0)(0)(0) << " (expect 0)" << std::endl;
  }
  

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
