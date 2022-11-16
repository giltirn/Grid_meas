#include<common/propagator_invert.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/lanczos.h>

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

  GridParallelRNG pRNG5rb(GridsD.FrbGrid);  
  pRNG5rb.SeedFixedIntegers(seeds4);

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
  
  LanczosParameters params;
  
  int Nstop = 30;
  int Nk = 30;
  int Np = 5;
  int ord = 80;
  RealD lo = 1.5;
  RealD hi = 88.0;
  params.alpha = sqrt(hi);
  params.beta = sqrt(lo);
  params.n_stop = Nstop;
  params.n_want = Nk;
  params.n_use = Nk+Np;
  params.ord = ord;
  params.tolerance = 1e-8;

  std::vector<RealD> eval_GP;
  std::vector<FermionFieldD> evec_GP;

  computeEigenvalues(eval_GP, evec_GP, params, GridsD.FGrid, GridsD.FrbGrid, Ud, *actions.action_d, pRNG5);


  std::vector<RealD> eval_Xconj;
  std::vector<FermionField1fD> evec_Xconj;

  computeEigenvalues(eval_Xconj, evec_Xconj, params, GridsD.FGrid, GridsD.FrbGrid, Ud, *actions.xconj_action_d, pRNG5);

  assert(eval_GP.size() == eval_Xconj.size());

  int Nevec = eval_GP.size();

  std::cout << "Comparing evals: " << std::endl;
  for(int i=0;i<Nevec;i++){
    std::cout << eval_Xconj[i] << " " << eval_GP[i] << " diff: " << eval_Xconj[i] - eval_GP[i] << std::endl;
    assert( fabs(eval_Xconj[i] - eval_GP[i]) < 1e-7 );
  }

  std::vector<FermionFieldD> evec_GP_rot(Nevec, GridsD.FrbGrid);
  
  std::cout << "Comparing eigenvectors: " << std::endl;
  for(int i=0;i<Nevec;i++){
    //We need to phase-rotate the regular evecs into X-conjugate vectors
    FermionField1fD v0 = PeekIndex<GparityFlavourIndex>(evec_GP[i],0);
    FermionField1fD v1 = PeekIndex<GparityFlavourIndex>(evec_GP[i],1);
    FermionField1fD Xv1star = Xmatrix()*conjugate(v1);
    ComplexD z = innerProduct(v0, Xv1star);
    ComplexD alpha = ComplexD(0.5)/z;

    FermionFieldD w = 1./sqrt(alpha) * evec_GP[i];
    evec_GP_rot[i] = w;

    std::cout << "Phase-rotated regular evec norm: " << norm2(w) << " vs original norm " << norm2(evec_GP[i]) << std::endl;

    //w should be X-conjugate; check
    FermionField1fD w0 = PeekIndex<GparityFlavourIndex>(w,0);
    FermionField1fD w1 = PeekIndex<GparityFlavourIndex>(w,1);
    FermionField1fD tmp = w1 + Xmatrix()*conjugate(w0);
    std::cout << "Check phase-rotated regular evec is X-conjugate (expect 0): " << norm2(tmp) << std::endl;
    assert(norm2(tmp) < 1e-8);

    tmp = w0 - evec_Xconj[i];
    FermionField1fD tmp2 = w0 + evec_Xconj[i];
    RealD nrmdiff = norm2(tmp);
    RealD nrmsum = norm2(tmp2);
    std::cout << "Evec " << i << " check X-conj Dop evec vs phase-rotated regular evec (expect diff or sum to be 0): difference: " << nrmdiff << " sum: " << nrmsum << std::endl;
    assert( nrmsum < 1e-6 || nrmdiff < 1e-6 );
  }

  std::cout << "Testing guessers" << std::endl;


  std::vector<FermionField1fF> evec_Xconj_f(Nevec, GridsF.FrbGrid);
  for(int i=0;i<Nevec;i++) precisionChange(evec_Xconj_f[i], evec_Xconj[i]);

  std::vector<FermionFieldF> evec_GP_rot_f(Nevec, GridsF.FrbGrid);
  for(int i=0;i<Nevec;i++) precisionChange(evec_GP_rot_f[i], evec_GP_rot[i]);

  std::vector<FermionFieldF> evec_GP_f(Nevec, GridsF.FrbGrid);
  for(int i=0;i<Nevec;i++) precisionChange(evec_GP_f[i], evec_GP[i]);

  {
    std::cout << "Arbitrary 2f source" << std::endl;   
    FermionFieldD src(GridsD.FrbGrid), tmp(GridsD.FrbGrid);
    gaussian(pRNG5rb, src);

    DeflatedGuesser<FermionFieldD> gGP(evec_GP, eval_GP);
    FermionFieldD guess_GP(GridsD.FrbGrid); guess_GP = Zero();
    gGP(src, guess_GP);

    RealD nrm2_GP = norm2(guess_GP);

    MixedPrecDeflatedGuesser<FermionFieldF,FermionFieldD> gGPmx(evec_GP_f, eval_GP);
    FermionFieldD guess_GPmx(GridsD.FrbGrid); guess_GPmx = Zero();
    gGPmx(src, guess_GPmx);
    
    tmp = guess_GPmx - guess_GP;
    std::cout << "GP mixed vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;
    
    DeflatedGuesser<FermionFieldD> gGP_rot(evec_GP_rot, eval_GP);
    FermionFieldD guess_GP_rot(GridsD.FrbGrid); guess_GP_rot = Zero();
    gGP_rot(src, guess_GP_rot);
    
    tmp = guess_GP_rot - guess_GP;
    std::cout << "GP rotated double vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;
    
    XconjDeflatedGuesser2fConvert gXconj(evec_Xconj, eval_Xconj);
    FermionFieldD guess_Xconj(GridsD.FrbGrid); guess_Xconj = Zero();
    gXconj(src, guess_Xconj);

    tmp = guess_Xconj - guess_GP;
    std::cout << "Xconj double (conv 2f) vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;

    XconjMixedPrecDeflatedGuesser2fConvert gXconjmx(evec_Xconj_f,eval_Xconj);
    FermionFieldD guess_Xconjmx(GridsD.FrbGrid); guess_Xconjmx = Zero();
    gXconjmx(src, guess_Xconjmx);

    tmp = guess_Xconjmx - guess_GP;
    std::cout << "Xconj mixed (conv 2f) vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;
  }

  {
    std::cout << "X-conjugate 2f source" << std::endl;   
    FermionField1fD src_1f(GridsD.FrbGrid);
    gaussian(pRNG5rb, src_1f);

    FermionFieldD src_2f(GridsD.FrbGrid), tmp(GridsD.FrbGrid);
    get2fXconjVector(src_2f, src_1f);

    DeflatedGuesser<FermionFieldD> gGP(evec_GP, eval_GP);
    FermionFieldD guess_GP(GridsD.FrbGrid); guess_GP = Zero();
    gGP(src_2f, guess_GP);

    RealD nrm2_GP = norm2(guess_GP);
    
    XconjDeflatedGuesser<FermionField1fD> gXconj(evec_Xconj, eval_Xconj);
    FermionField1fD guess_Xconj(GridsD.FrbGrid); guess_Xconj = Zero();
    gXconj(src_1f, guess_Xconj);
    FermionFieldD guess_Xconj_2f(GridsD.FrbGrid); 
    get2fXconjVector(guess_Xconj_2f, guess_Xconj);

    tmp = guess_Xconj_2f - guess_GP;
    std::cout << "Xconj double (native) vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;

    XconjMixedPrecDeflatedGuesser<FermionField1fF,FermionField1fD> gXconjmx(evec_Xconj_f, eval_Xconj);
    FermionField1fD guess_Xconjmx(GridsD.FrbGrid); guess_Xconjmx = Zero();
    gXconjmx(src_1f, guess_Xconjmx);
    FermionFieldD guess_Xconjmx_2f(GridsD.FrbGrid); 
    get2fXconjVector(guess_Xconjmx_2f, guess_Xconjmx);

    tmp = guess_Xconjmx_2f - guess_GP;
    std::cout << "Xconj mixed (native) vs GP double (expect 0): " << norm2(tmp)/nrm2_GP << std::endl;
  }
    
  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
