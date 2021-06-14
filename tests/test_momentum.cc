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

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(UGridD); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  //Gparity in 1 direction
  assert(Nd == 4);
  std::vector<int> dirs4(4, 0);
  dirs4[1] = 1;
  
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = 0;
  Params.twists[1] = 1; //GPBC in Y
  Params.twists[Nd-1] = 1; //APBC in time direction

  {
    //Test of physical momentum
    std::vector<int> pbase = {0,1,1}; 
    std::vector<double> pphys = getPhysicalMomentum(pbase);
    assert(pphys.size() == 3);
    assert(pphys[0] == 0.);
    assert(pphys[1] == M_PI/2./latt[1]);
    assert(pphys[2] == 2*M_PI/latt[2]);
    std::cout << GridLogMessage << "Passed physical momentum test" << std::endl;
  }
  

  {
    //Test phase at site
    Coordinate x(std::vector<int>({1,1,1,1}));
    std::vector<double> p = {0, M_PI/2., 0};
    ComplexD z = phase(p, x); //exp( -i pi/2 ) = -i
    assert( fabs(z.real()) < 1e-6 );
    assert( fabs(z.imag()+1.) < 1e-6 );
    std::cout << GridLogMessage << "Passed site phase test" << std::endl;
  }

  {
    //Test phase field
    std::vector<double> p = {0, M_PI/2., 0};    
    LatticeComplexD pf = phaseField(p, UGridD);

    Coordinate c1(std::vector<int>({1,1,1,1}));

    ComplexD vc1;
    peekSite(vc1, pf, c1); //expect -i
    
    assert( fabs(vc1.real()) < 1e-6 );
    assert( fabs(vc1.imag()+1.) < 1e-6 );

    Coordinate c2(std::vector<int>({0,2,0,0}));

    ComplexD vc2;
    peekSite(vc2, pf, c2); //expect -1

    assert( fabs(vc2.real()+1.) < 1e-6 );
    assert( fabs(vc2.imag()) < 1e-6 );
    
    std::cout << GridLogMessage << "Passed phase field test" << std::endl;
  }

  {
    //Test projectors
    //psi_+  :   p = (..., -7, -3, 1, 5, 9, ...) pi/2L
    //psi_-  :   p = (..., -9, -5, -1, 3, 7, ...) pi/2L

    std::vector<int> p(3,0);

    p[1] = 9;
    GparityFlavour::Algebra alg = getProjector(p);
    assert(alg == GparityFlavour::Algebra::ProjPlus);
    
    p[1] = -7;
    alg = getProjector(p);
    assert(alg == GparityFlavour::Algebra::ProjPlus);


    p[1] = 7;
    alg = getProjector(p);
    assert(alg == GparityFlavour::Algebra::ProjMinus);
  
    p[1] = -9;
    alg = getProjector(p);
    assert(alg == GparityFlavour::Algebra::ProjMinus);
    
    std::cout << GridLogMessage << "Passed projector test" << std::endl;
  }    

  {
    //Test momstr
    std::vector<int> p = {1,2,3};
    std::string pstr = momstr(p);
    assert(pstr == "123");
    p = {-1,2,-3};
    pstr = momstr(p);
    assert(pstr == "_12_3");

    std::cout << GridLogMessage << "Passed momstr test" << std::endl;
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
