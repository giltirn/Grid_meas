#include<common/sources.h>


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

  {
    std::cout << GridLogMessage << "Testing point source" << std::endl;

    SCFmatrixD one = unitSiteSrc();
    assert(norm2(one) == RealD(Ngp*Ns*Nc) );
    assert( TensorRemove(trace(one)) == ComplexD(Ngp*Ns*Nc) );

    SCFmatrixD zero = Zero();
    Coordinate coord( std::vector<int>({1,0,1,0}) );

    LatticeSCFmatrixD pnt = pointSource(coord, UGridD);

    SCFmatrixD v;
    peekSite(v, pnt, coord);
    
    std::cout << GridLogMessage << "At " << coord << " got " << v << " expect " << one << std::endl;
    assert(v == one);

    //Choose 20 random points and check
    std::default_random_engine generator;
    Coordinate c(Nd);

    for(int i=0; i<20; i++){

      randomCoordinate(c, latt, 4, generator, UGridD);
      peekSite(v, pnt, c);
      std::cout << GridLogMessage << "At " << c << " got " << v << " expect " << ( c==coord ? one : zero ) << std::endl;

      if(c == coord) assert(v == one);
      else assert(v == zero);
    }
    
    std::cout << GridLogMessage << "Point test passed" << std::endl;
  }


  {
    std::cout << GridLogMessage << "Testing Z2 source" << std::endl;

    SCFmatrixD one = unitSiteSrc();
    SCFmatrixD mone = -one;
    SCFmatrixD zero = Zero();
    Coordinate coord( std::vector<int>({1,0,1,0}) );

    LatticeSCFmatrixD wall =Z2wallSource(1, pRNG, UGridD);

    SCFmatrixD v;

    //Choose 20 random points on the wall and check
    std::default_random_engine generator;
    Coordinate c(Nd);
    int tsrc=1;
    c[3] = tsrc;
    
    for(int i=0; i<20; i++){
      randomCoordinate(c, latt, 3, generator, UGridD);
      peekSite(v, wall, c);

      int mdim = Ngp*Ns*Nc;
      assert( norm2(v) == mdim );
      ComplexD tr = TensorRemove(trace(v));
      std::cout << GridLogMessage << "At " << c << " got trace " << tr << std::endl;
      tr = tr / RealD(mdim);
      assert( tr == ComplexD(1.0) || tr == ComplexD(-1.0) );
      int sign = tr == ComplexD(1.0) ? 1 : -1;
      
      std::cout << " and sign " << sign << std::endl;
    }

    //20 points off the wall
    for(int i=0; i<20; i++){
      randomCoordinate(c, latt, 4, generator, UGridD);
      if(c[3] == tsrc) 
	c[3] == (tsrc + 1) % latt[3];

      peekSite(v, wall, c);

      ComplexD tr = TensorRemove(trace(v));
      std::cout << GridLogMessage << "At " << c << " got trace " << tr << " expect 0" << std::endl;

      assert(v == zero);
    }

    std::cout << GridLogMessage << "Z2Wall test passed" << std::endl;
  }

  {
    std::cout << GridLogMessage << "Testing wall source" << std::endl;

    SCFmatrixD one = unitSiteSrc();
    SCFmatrixD zero = Zero();
    Coordinate coord( std::vector<int>({1,0,1,0}) );

    LatticeSCFmatrixD wall =wallSource(1, UGridD);

    SCFmatrixD v;

    //Choose 20 random points on the wall and check
    std::default_random_engine generator;
    Coordinate c(Nd);
    int tsrc=1;
    c[3] = tsrc;
    
    for(int i=0; i<20; i++){
      randomCoordinate(c, latt, 3, generator, UGridD);
      peekSite(v, wall, c);

      int mdim = Ngp*Ns*Nc;
      assert( norm2(v) == mdim );
      ComplexD tr = TensorRemove(trace(v));
      std::cout << GridLogMessage << "At " << c << " got trace " << tr << std::endl;
      tr = tr / RealD(mdim);
      assert( tr == ComplexD(1.0) );
    }

    //20 points off the wall
    for(int i=0; i<20; i++){
      randomCoordinate(c, latt, 4, generator, UGridD);
      if(c[3] == tsrc) 
	c[3] == (tsrc + 1) % latt[3];

      peekSite(v, wall, c);

      ComplexD tr = TensorRemove(trace(v));
      std::cout << GridLogMessage << "At " << c << " got trace " << tr << " expect 0" << std::endl;

      assert(v == zero);
    }

    std::cout << GridLogMessage << "Wall test passed" << std::endl;
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
