#include<common/propagator_invert.h>
#include<common/utils.h>


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

  LatticeSCFmatrixD rsrc(UGridD);
  gaussian(pRNG, rsrc);

  std::default_random_engine generator;
  {
    std::cout << GridLogMessage << "Checking column peek/poke" << std::endl;
    LatticeSCFmatrixD recon(UGridD);
    recon = Zero();

    Coordinate coord(Nd);
    typename GridTypeMapper<GparityWilsonImplD::SiteSpinor>::scalar_object ss;    
    SCFmatrixD full_site;
    FermionFieldD tmp(UGridD);
    for(int f=0;f<Ngp;f++){
      for(int s=0;s<Ns;s++){
	for(int c=0;c<Nc;c++){
	  tmp = extractColumn(rsrc, f, s, c);
	  
	  //Try 10 different random sites
	  for(int i=0;i<10;i++){
	    randomCoordinate(coord, latt, 4, generator, UGridD);
	    peekSite(full_site, rsrc, coord);
	    peekSite(ss, tmp, coord);

	    for(int rf=0;rf<Ngp;rf++){
	      for(int rs=0;rs<Ns;rs++){
		for(int rc=0;rc<Nc;rc++){		  
		  assert( ss(rf)(rs)(rc)== full_site(rf,f)(rs,s)(rc,c) );
		}
	      }
	    }
	    std::cout << "Passed column check " <<  f<< " " <<  s<< " " << c << " on site " << coord << std::endl;
	  }
	  insertColumn(recon, tmp, f,s,c);
	}
      }
    }
    LatticeSCFmatrixD diff = recon - rsrc;
    RealD diff_nrm2= norm2(diff);
    std::cout << GridLogMessage << "Diff between original and reconstructed src " << diff_nrm2 << std::endl;

    assert( fabs(diff_nrm2) < 1e-10 );
    std::cout << GridLogMessage << "Column peek/poke test passed" << std::endl;
  }
  
  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}

	  

