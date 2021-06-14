#include<common/array_io.h>

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
    //Test binary io of real vector
    std::vector<RealD> in(30);
    for(int i=0;i<30;i++) in[i] = M_PI * (i+1);
    binaryWriteArray("test.bin", in);

    std::vector<RealD> out;
    binaryReadArray(out, "test.bin");
    
    for(int i=0;i<30;i++)
      assert(out[i] == in[i]);

    std::cout << GridLogMessage << "Passed binary io of real vector test" << std::endl;
  }
  
  {
    //Test ASCII write
    std::vector<RealD> in(30);
    for(int i=0;i<30;i++) in[i] = M_PI * (i+1);

    asciiWriteArray(in, "test", 0);
    assert(fileExists("test.0"));

    std::ifstream f("test.0");
    std::vector<RealD> out(30);
    for(int i=0;i<30;i++)
      f >> out[i];

    for(int i=0;i<30;i++)
      assert(fabs(out[i] - in[i]) < 1e-10 );

    std::cout << GridLogMessage << "Passed ascii write of real vector test" << std::endl;
  }

  {
    //Test array of fields
    std::vector<DomainWallFermionD::FermionField> fields(10, FGridD);
    for(int i=0;i<10;i++)
      gaussian(pRNG, fields[i]);

    writeFieldArray("test.bin", fields);
    
    std::vector<DomainWallFermionD::FermionField> rdfields(10, FGridD);
    readFieldArray(rdfields, "test.bin");

    DomainWallFermionD::FermionField diff(FGridD);
    for(int i=0;i<10;i++){
      diff = rdfields[i] - fields[i];
      RealD n = norm2(diff);
      assert(n < 1e-12);
    }
    std::cout << GridLogMessage << "Passed binary write of field vector test" << std::endl;
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}


