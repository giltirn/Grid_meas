#include<common/guesser.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/lanczos.h>
#include<common/field_utils.h>

using namespace Grid;
using namespace GridMeas;

void test(const FermionFieldD &sol2f, const FermionFieldD &sol1f, const std::string &descr, const double tol){
  RealD sz1 = norm2(sol1f);
  RealD sz2 = norm2(sol2f);
  FermionFieldD tmp = sol1f - sol2f;
  RealD diff = norm2(tmp);
  RealD reldiff = 2*diff/(sz1+sz2);
  
  std::cout << descr << ": " << sz1 << " " << sz2 << " diff " << diff << " reldiff " << reldiff << " (expect 0)" << std::endl;
  assert(reldiff < tol);
}
void test(const FermionFieldD &sol2f, const FermionField1fD &sol1f, const std::string &descr, const double tol){
  FermionFieldD sol1f_conv(sol1f.Grid());
  get2fXconjVector(sol1f_conv, sol1f);
  test(sol2f, sol1f_conv, descr, tol);
}


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

  //Generate some fake X-conjugate eigenvectors and eigenvalues
  int Nevec = 20;
  std::vector<FermionField1fD> evec_1f(Nevec, GridsD.FrbGrid);
  std::vector<FermionFieldD> evec_2f(Nevec, GridsD.FrbGrid);
  std::vector<RealD> eval(Nevec);
  for(int i=0;i<Nevec;i++){
    gaussian(pRNG,evec_1f[i]);
    get2fXconjVector(evec_2f[i], evec_1f[i]);
    gaussian(sRNG, eval[i]);
    eval[i] = 0.1 + fabs(eval[i]); //positive-definite
  }

  ////////////////////////////////////////////////////////////////////////////
  //       Double precision
  ////////////////////////////////////////////////////////////////////////////

  std::cout << "Double precision" << std::endl;

  //Create the guessers
  XconjDeflatedGuesser<FermionField1fD> guesser_1f(evec_1f,eval);
  DeflatedGuesser<FermionFieldD> guesser_2f(evec_2f,eval);
  XconjDeflatedGuesser2fConvert guesser_2f_econv(evec_1f,eval);
  
  //Create sources
  FermionField1fD src_1f(GridsD.FrbGrid);
  gaussian(pRNG,src_1f);

  FermionFieldD src_2f(GridsD.FrbGrid);
  get2fXconjVector(src_2f, src_1f);
  std::cout << "Source norm " << norm2(src_2f) << std::endl;

  FermionFieldD src_2f_nonXconj(GridsD.FrbGrid);
  gaussian(pRNG,src_2f_nonXconj);
  
  //Test
  FermionFieldD sol_2f(GridsD.FrbGrid); sol_2f = Zero();
  guesser_2f(src_2f, sol_2f);

  FermionField1fD sol_1f(GridsD.FrbGrid); sol_1f = Zero();
  guesser_1f(src_1f, sol_1f);

  test(sol_2f, sol_1f, "2f double vs 1f double (native) with Xconj src", 1e-5); 

  FermionFieldD sol_2f_econv(GridsD.FrbGrid); sol_2f_econv = Zero();
  guesser_2f_econv(src_2f, sol_2f_econv);
  test(sol_2f, sol_2f_econv, "2f double vs 1f double (conv) with Xconj src", 1e-5); 
  
  FermionFieldD sol_2f_nonXconj(GridsD.FrbGrid); sol_2f_nonXconj = Zero();
  guesser_2f(src_2f_nonXconj, sol_2f_nonXconj);
  
  sol_2f_econv = Zero();
  guesser_2f_econv(src_2f_nonXconj, sol_2f_econv);
  test(sol_2f_nonXconj, sol_2f_econv, "2f double vs 1f double (conv) with non-Xconj src", 1e-5); 

  ////////////////////////////////////////////////////////////////////////////
  //       Single precision
  ////////////////////////////////////////////////////////////////////////////

  //Create converted evecs
  std::vector<FermionField1fF> evec_1f_f(Nevec, GridsF.FrbGrid);
  std::vector<FermionFieldF> evec_2f_f(Nevec, GridsF.FrbGrid);
  for(int i=0;i<Nevec;i++){
    precisionChange(evec_1f_f[i],evec_1f[i]);
    precisionChange(evec_2f_f[i],evec_2f[i]);
  }

  XconjMixedPrecDeflatedGuesser<FermionField1fF, FermionField1fD> guesser_mx_1f(evec_1f_f,eval);
  MixedPrecDeflatedGuesser<FermionFieldF, FermionFieldD> guesser_mx_2f(evec_2f_f,eval);
  XconjMixedPrecDeflatedGuesser2fConvert guesser_mx_2f_econv(evec_1f_f,eval);
    
  FermionFieldD sol_2f_mx(GridsD.FrbGrid); sol_2f_mx = Zero();
  guesser_mx_2f(src_2f, sol_2f_mx);
  test(sol_2f, sol_2f_mx, "2f double vs 2f mixed with Xconj src", 1e-5); 

  sol_1f = Zero();
  guesser_mx_1f(src_1f, sol_1f);
  test(sol_2f, sol_1f, "2f double vs 1f mixed (native) with Xconj src", 1e-5); 

  sol_2f_econv = Zero();
  guesser_mx_2f_econv(src_2f, sol_2f_econv);
  test(sol_2f, sol_2f_econv, "2f double vs 1f mixed (conv) with Xconj src", 1e-5); 

  sol_2f_econv = Zero();
  guesser_mx_2f_econv(src_2f_nonXconj, sol_2f_econv);
  test(sol_2f_nonXconj, sol_2f_econv, "2f double vs 1f mixed (conv) with non-Xconj src", 1e-5); 


  std::cout << "Checking getGuesser" << std::endl;
  {
    LinearFunction<FermionFieldD>* gptr = getGuesser<FermionFieldD,FermionFieldD>(eval, evec_2f);
    DeflatedGuesser<FermionFieldD>* test = dynamic_cast<DeflatedGuesser<FermionFieldD>* >(gptr);
    assert(test != nullptr);
  }
  {
    LinearFunction<FermionFieldD>* gptr = getGuesser<FermionFieldD,FermionFieldF>(eval, evec_2f_f);
    MixedPrecDeflatedGuesser<FermionFieldF,FermionFieldD>* test = dynamic_cast<MixedPrecDeflatedGuesser<FermionFieldF,FermionFieldD>* >(gptr);
    assert(test != nullptr);
  }
  {
    LinearFunction<FermionField1fD>* gptr = getGuesser<FermionField1fD,FermionField1fD>(eval, evec_1f);
    XconjDeflatedGuesser<FermionField1fD>* test = dynamic_cast<XconjDeflatedGuesser<FermionField1fD>* >(gptr);
    assert(test != nullptr);
  }
  {
    LinearFunction<FermionField1fD>* gptr = getGuesser<FermionField1fD,FermionField1fF>(eval, evec_1f_f);
    XconjMixedPrecDeflatedGuesser<FermionField1fF,FermionField1fD>* test = dynamic_cast<XconjMixedPrecDeflatedGuesser<FermionField1fF,FermionField1fD>* >(gptr);
    assert(test != nullptr);
  }
  {
    LinearFunction<FermionFieldD>* gptr = getGuesser<FermionFieldD,FermionField1fD>(eval, evec_1f);
    XconjDeflatedGuesser2fConvert* test = dynamic_cast<XconjDeflatedGuesser2fConvert* >(gptr);
    assert(test != nullptr);
  }
  {
    LinearFunction<FermionFieldD>* gptr = getGuesser<FermionFieldD,FermionField1fF>(eval, evec_1f_f);
    XconjMixedPrecDeflatedGuesser2fConvert* test = dynamic_cast<XconjMixedPrecDeflatedGuesser2fConvert* >(gptr);
    assert(test != nullptr);
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
