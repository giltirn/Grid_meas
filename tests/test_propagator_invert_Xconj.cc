#include<common/propagator_invert_Xconj.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/sources.h>

using namespace Grid;
using namespace GridMeas;

  
  //Invert using the X-conjugate action for flavor diagonal sources
  //NOTE: the source *must* have flavor structure  
  //|  A       0   |
  //|  0   -XA*X |
  //which includes
  // | phi   0   |
  // |  0   phi* |  
  //where phi is a spin-color matrix that commutes with X=C g5  
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertGenXconjDiag(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
			       const LatticeSCFmatrixD &msrc, FermionActionD &xconj_action_d, FermionActionF &xconj_action_f, 
			       const MixedCGargs &args,
			       bool do_midprop,
			       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){

    LatticePropagatorD A = PeekIndex<GparityFlavourIndex>(msrc,0,0);
    LatticePropagatorD MXA_sol(msrc.Grid()), MXA_mid(msrc.Grid());
    mixedPrecInvertGen(MXA_sol, MXA_mid, A, xconj_action_d, xconj_action_f, args, do_midprop, evals, evecs);

    std::vector<LatticePropagatorD> prop_fcol(2, msrc.Grid());
    prop_fcol[0] = mulPplusRight(MXA_sol);
    prop_fcol[1] = mulPminusRight(MXA_sol);

    get2fXconjMatrix(prop, prop_fcol[0],prop_fcol[1]);
    prop = mulHdagRight(prop);

    if(do_midprop){
      prop_fcol[0] = mulPplusRight(MXA_mid);
      prop_fcol[1] = mulPminusRight(MXA_mid);

      get2fXconjMatrix(midprop, prop_fcol[0],prop_fcol[1]);
      midprop = mulHdagRight(midprop);
    }
  }


  //No midprop, with evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  LatticeSCFmatrixD mixedPrecInvertXconjDiag(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, const MixedCGargs &args,
					 std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    LatticeSCFmatrixD tmp(msrc.Grid()), prop(msrc.Grid());
    mixedPrecInvertGenXconjDiag(prop,tmp,msrc,action_d,action_f,args,false,evals,evecs);
    return prop;
  }
  //No midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  LatticeSCFmatrixD mixedPrecInvertXconjDiag(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, const MixedCGargs &args){
    return mixedPrecInvertXconjDiag(msrc,action_d,action_f,args,(std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }


  //With midprop and evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertWithMidPropXconjDiag(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       const MixedCGargs &args,
				       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    mixedPrecInvertGenXconjDiag(prop,midprop,msrc,action_d,action_f,args,true,evals,evecs);
  }    
  //With midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertWithMidPropXconjDiag(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       const MixedCGargs &args){
    mixedPrecInvertWithMidPropXconjDiag(prop, midprop, msrc, action_d, action_f, args, (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }
 



void projectXconj(FermionFieldD &out, const FermionFieldD &f){
  FermionField1fD phi1 = PeekIndex<GparityFlavourIndex>(f,0);
  FermionField1fD phi2 = PeekIndex<GparityFlavourIndex>(f,1);
  FermionField1fD rho = 0.5*( phi1 + Xmatrix()*conjugate(phi2));
  get2fXconjVector(out,rho);
}
void projectXbarConj(FermionFieldD &out, const FermionFieldD &f){
  FermionField1fD phi1 = PeekIndex<GparityFlavourIndex>(f,0);
  FermionField1fD phi2 = PeekIndex<GparityFlavourIndex>(f,1);
  FermionField1fD rho = 0.5*( phi1 - Xmatrix()*conjugate(phi2));
  get2fXbarConjVector(out,rho);
}

//FermionFieldD 
template<typename FieldType>
void multXiLeft(FieldType &to, const FieldType &from){
  static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
  to = ComplexD(0,-1)*(sigma2*(Xmatrix()*from));
}
template<typename FieldType>
void multXiRight(FieldType &to, const FieldType &from){
  static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
  to = ComplexD(0,-1)*((from*sigma2)*Xmatrix());
}



// template<typename FermionActionD, typename FermionActionF>
// void invertFullyDilutedSourceXconj(LatticeSCFmatrixD &prop,
// 				   const LatticeComplexD &rho,
// 				   FermionActionD &xconj_action_d, FermionActionF &xconj_action_f,
// 				   const MixedCGargs &args
// 				   ){
//   SchurDiagMooeeOperator<FermionActionD,FermionField1fD> hermop_d(xconj_action_d);
//   SchurDiagMooeeOperator<FermionActionF,FermionField1fF> hermop_f(xconj_action_f);
  
//   LinearFunction<FermionField1fD>* mcg = MixedCGfactory(args, xconj_action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
//   LinearFunctionWrapper<FermionField1fD> mcg_wrap(*mcg);

//   LinearFunction<FermionField1fD> *guesser = nullptr;
//   SchurRedBlackDiagMooeeSolve<FermionField1fD> solver(mcg_wrap);
  
//   GridBase* FGridD = xconj_action_d.FermionGrid();

//   conformable(rho.Grid(), prop.Grid());

//   FermionField1fD src_1f_5d(FGridD), sol_1f_5d(FGridD);
//   FermionField1fD src_1f_4d(rho.Grid()), sol_1f_4d_X(rho.Grid()), sol_1f_4d_Xbar(rho.Grid());
//   FermionFieldD sol_2f_4d_X(rho.Grid()), sol_2f_4d_Xbar(rho.Grid()), tmp_2f_1(rho.Grid()), tmp_2f_2(rho.Grid()), tmp_2f_3(rho.Grid());

//   LatticeColourVectorD tmp_c(rho.Grid());
//   LatticeSCFmatrixD tmp_scfmat(rho.Grid());
//   tmp_scfmat = Zero();
  
//   for(int s=0;s<Ns;s++){
//     for(int c=0;c<Nc;c++){
//       tmp_c = Zero();
//       PokeIndex<ColourIndex>(tmp_c,rho,c);
//       src_1f_4d = Zero();
//       PokeIndex<SpinIndex>(src_1f_4d, tmp_c, s);
      
//       //sol_1f_4d_X = M^-1_X rho
//       xconj_action_d.ImportPhysicalFermionSource(src_1f_4d, src_1f_5d);
//       solver(xconj_action_d, src_1f_5d, sol_1f_5d);
//       xconj_action_d.ExportPhysicalFermionSolution(sol_1f_5d, sol_1f_4d_X);

//       //Also need   M^-1_Xbar rho = -i M^-1_X (i rho)
//       src_1f_4d = ComplexD(0,1)*src_1f_4d;
//       xconj_action_d.ImportPhysicalFermionSource(src_1f_4d, src_1f_5d);
//       solver(xconj_action_d, src_1f_5d, sol_1f_5d);
//       xconj_action_d.ExportPhysicalFermionSolution(sol_1f_5d, sol_1f_4d_Xbar);
//       sol_1f_4d_Xbar = ComplexD(0,-1)*sol_1f_4d_Xbar;

//       //Now we can form the f=0 column
//       get2fXconjVector(sol_2f_4d_X,sol_1f_4d_X);
//       get2fXbarConjVector(sol_2f_4d_Xbar,sol_1f_4d_Xbar);
      
//       tmp_2f_1 = sol_2f_4d_X+sol_2f_4d_Xbar;
//       insertColumn(prop, tmp_2f_1, 0,s,c);
      
//       //Store Xi sol_2f_4d_X^* - Xi sol_2f_4d_Xbar^* as the appropriate column of a temp matrix as we need to right-multiply by X
//       multXiLeft(tmp_2f_1, conjugate(sol_2f_4d_X));
//       multXiLeft(tmp_2f_2, -conjugate(sol_2f_4d_Xbar));
//       tmp_2f_3 = tmp_2f_1 + tmp_2f_2;
//       insertColumn(tmp_scfmat, tmp_2f_3, 1, s,c);
//     }
//   }
//   tmp_scfmat = tmp_scfmat * Xmatrix(); //manipulate spin columns
//   for(int s=0;s<Ns;s++){
//     for(int c=0;c<Nc;c++){
//       tmp_2f_1 = extractColumn(tmp_scfmat, 1,s,c);
//       insertColumn(prop, tmp_2f_1, 1,s,c);
//     }
//   }

// }




template<typename FermionActionD, typename FermionActionF>
void invertFullyDilutedSourceXconj(LatticeSCFmatrixD &prop,
				   const LatticeComplexD &rho,
				   FermionActionD &xconj_action_d, FermionActionF &xconj_action_f,
				   const MixedCGargs &args
				   ){
  //First solve
  //M^-1 | rho  0   |
  //     |  0  rho* |

  LatticeSCFmatrixD src_cconj(rho.Grid());
  
  {
    size_t Nsimd = FermionFieldD::vector_type::Nsimd();
    autoView( rho_v, rho, AcceleratorRead);
    autoView( src_cconj_v, src_cconj, AcceleratorWrite);
    accelerator_for( i, rho_v.size(), Nsimd, {
	auto site_rho = rho_v(i);
	auto site_into = src_cconj_v(i);
	for(int f=0;f<Ngp;f++)
	  for(int s=0;s<Ns;s++)
	    for(int c=0;c<Nc;c++)
	      site_into(f,f)(s,s)(c,c) = (f==0 ? site_rho()()() : conjugate(site_rho()()()) );
	coalescedWrite(src_cconj_v[i], site_into);
      });
  }
  LatticeSCFmatrixD psi = mixedPrecInvertXconj(src_cconj, xconj_action_d, xconj_action_f, args);
  
  static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
  
  LatticeSCFmatrixD &tmp = src_cconj; //reuse
  multXiLeft(tmp,conjugate(psi));
  LatticeSCFmatrixD Xi_psistar_Xi(rho.Grid());
  multXiRight(Xi_psistar_Xi,tmp);

  prop = psi + psi*sigma3 + Xi_psistar_Xi - Xi_psistar_Xi*sigma3;
  prop = ComplexD(0.5)*prop;
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

  Actions actions(ActionType::Mobius, Params, 0.01, 2.0, Ud, GridsD, Uf, GridsF);
  
  //Should work with any matrix of the form
  //|  A       B   |
  //| XB*X   -XA*X |

  LatticePropagator A(GridsD.UGrid);
  gaussian(pRNG, A);

  LatticePropagator B(GridsD.UGrid);
  gaussian(pRNG, B);
  
  LatticePropagator C = Xmatrix()*(conjugate(B)*Xmatrix());
  LatticePropagator D = -(Xmatrix()*(conjugate(A)*Xmatrix()));
  
  LatticeSCFmatrixD src(GridsD.UGrid);
  src = Zero();

  PokeIndex<GparityFlavourIndex>(src, A, 0,0);
  PokeIndex<GparityFlavourIndex>(src, B, 0,1);
  PokeIndex<GparityFlavourIndex>(src, C, 1,0);
  PokeIndex<GparityFlavourIndex>(src, D, 1,1);
 
  //The Xconj-diag implementation works only for matrices of the form
  //|  A       0   |
  //|  0   -XA*X   |

  LatticeSCFmatrixD src_diag(GridsD.UGrid);
  src_diag = Zero();

  PokeIndex<GparityFlavourIndex>(src_diag, A, 0,0);
  PokeIndex<GparityFlavourIndex>(src_diag, D, 1,1);

  //Test left and right multiplication by H
  {
    LatticeSCFmatrixD r1 = mulHdagLeft(mulHLeft(src));
    LatticeSCFmatrixD r2 = mulHLeft(mulHdagLeft(src));
    LatticeSCFmatrixD r3 = mulHdagRight(mulHRight(src));
    LatticeSCFmatrixD r4 = mulHRight(mulHdagRight(src));
    LatticeSCFmatrixD tmp(src.Grid());
    tmp = r1 - src;
    RealD diff1 = norm2(tmp);
    tmp = r2 - src;
    RealD diff2 = norm2(tmp);
    tmp = r3 - src;
    RealD diff3 = norm2(tmp);
    tmp = r4 - src;
    RealD diff4 = norm2(tmp);

    std::cout << "Test H and H^-1 left and right multiplication (expect 0): " << diff1 << " " << diff2 << " " << diff3 << " " << diff4 << std::endl;
    assert(diff1 < 1e-8);
    assert(diff2 < 1e-8);
    assert(diff3 < 1e-8);
    assert(diff4 < 1e-8);
  }
  //Test src*H has X-conjugate columns
  {
    std::cout << "Test src*H has X-conjugate columns:" << std::endl;
    LatticeSCFmatrixD r = mulHRight(src);
    typedef columnOps<FermionFieldD> cop;
    for(int c=0;c<24;c++){
      FermionFieldD col = cop::extractColumn(r, c);
      FermionField1fD col_0 = PeekIndex<GparityFlavourIndex>(col,0);
      FermionField1fD col_1 = PeekIndex<GparityFlavourIndex>(col,1);
      FermionField1fD tmp = col_1 + Xmatrix()*conjugate(col_0);
      std::cout << c << " (expect 0): " << norm2(tmp) << std::endl;
      assert(norm2(tmp) < 1e-10);
    }
  }

  MixedCGargs cg_args;
  cg_args.algorithm = MixedCGalgorithm::ReliableUpdateCG;

#if 0
  LatticeSCFmatrixD sol_GP(GridsD.UGrid), midsol_GP(GridsD.UGrid);
  mixedPrecInvertGen(sol_GP, midsol_GP, src, *actions.action_d, *actions.action_f, cg_args, true, 
		     (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);

  LatticeSCFmatrixD sol_Xconj(GridsD.UGrid), midsol_Xconj(GridsD.UGrid);
  mixedPrecInvertGenXconj(sol_Xconj, midsol_Xconj, src, *actions.xconj_action_d, *actions.xconj_action_f, cg_args, true, 
  			  (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);

  LatticeSCFmatrixD sol_GP_diag(GridsD.UGrid), midsol_GP_diag(GridsD.UGrid);
  mixedPrecInvertGen(sol_GP_diag, midsol_GP_diag, src_diag, *actions.action_d, *actions.action_f, cg_args, true, 
		     (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);

  LatticeSCFmatrixD sol_Xconj_diag(GridsD.UGrid), midsol_Xconj_diag(GridsD.UGrid);
  mixedPrecInvertGenXconjDiag(sol_Xconj_diag, midsol_Xconj_diag, src_diag, *actions.xconj_action_d, *actions.xconj_action_f, cg_args, true, 
			      (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);

  LatticeSCFmatrixD diff = sol_Xconj - sol_GP;
  std::cout << "Prop: Norm of difference GP vs Xconj (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);
  diff = midsol_Xconj - midsol_GP;
  std::cout << "Midprop: Norm of difference GP vs Xconj (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);

  diff = sol_Xconj_diag - sol_GP_diag;
  std::cout << "Prop: Norm of difference GP vs Xconj for diagonal source special case (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);
  diff = midsol_Xconj_diag - midsol_GP_diag;
  std::cout << "Midprop: Norm of difference GP vs Xconj  for diagonal source special case (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);
#endif



  //Test the strategy for computing sources of the form
  //| rho   0  |
  //|  0   rho |
  
  LatticeSCFmatrixD prop_test(GridsD.UGrid);
  LatticeComplexD rho(GridsD.UGrid);
  gaussian(pRNG, rho);
  
  {
    //First test the complex conjugate relation
    
    LatticeSCFmatrixD src_1(rho.Grid());
    
    {
      size_t Nsimd = FermionFieldD::vector_type::Nsimd();
      autoView( rho_v, rho, AcceleratorRead);
      autoView( src_1, src_1, AcceleratorWrite);
      accelerator_for( i, rho_v.size(), Nsimd, {
	  auto site_rho = rho_v(i);
	  auto site_into = src_1_v(i);
	  for(int f=0;f<Ngp;f++)
	    for(int s=0;s<Ns;s++)
	      for(int c=0;c<Nc;c++)
		site_into(f,f)(s,s)(c,c) = (f==0 ? site_rho()()() : conjugate(site_rho()()()) );
	  coalescedWrite(src_1_v[i], site_into);
	});
    }
    LatticeSCFmatrixD src_2(rho.Grid());
    {
      size_t Nsimd = FermionFieldD::vector_type::Nsimd();
      autoView( rho_v, rho, AcceleratorRead);
      autoView( src_2, src_2, AcceleratorWrite);
      accelerator_for( i, rho_v.size(), Nsimd, {
	  auto site_rho = rho_v(i);
	  auto site_into = src_2_v(i);
	  for(int s=0;s<Ns;s++)
	    for(int c=0;c<Nc;c++){
	      site_into(0,1)(s,s)(c,c) = -conjugate(site_rho);
	      site_into(1,0)(s,s)(c,c) = conjugate(site_rho);
	    }
	  coalescedWrite(src_2_v[i], site_into);
	});
    }
    src_2 = src_2 * Xmatrix();

    LatticeSCFmatrixD psi_1 = mixedPrecInvert(src_1, *actions.action_d, *actions.action_f, args);
    LatticeSCFmatrixD psi_2 = mixedPrecInvert(src_2, *actions.action_d, *actions.action_f, args);
    LatticeSCFmatrixD tmp(rho.grid());
    multXiLeft(tmp, conjugate(psi_1));
    
    



  }



#if 0
  invertFullyDilutedSourceXconj(prop_test,rho,*actions.xconj_action_d, *actions.xconj_action_f,cg_args);
  
  LatticeSCFmatrixD src_test(GridsD.UGrid);
  src_test = Zero();
  {
    size_t Nsimd = FermionFieldD::vector_type::Nsimd();
    autoView( rho_v, rho, AcceleratorRead);
    autoView( src_test_v, src_test, AcceleratorWrite);
    accelerator_for( i, rho_v.size(), Nsimd, {
	auto site_rho = rho_v(i);
	auto site_into = src_test_v(i);
	for(int f=0;f<Ngp;f++)
	  for(int s=0;s<Ns;s++)
	    for(int c=0;c<Nc;c++)
	      site_into(f,f)(s,s)(c,c) = site_rho()()();
	coalescedWrite(src_test_v[i], site_into);
      });
  }
  
  LatticeSCFmatrixD prop_ref = mixedPrecInvert(src_test, *actions.action_d, *actions.action_f, cg_args);
  diff = prop_test - prop_ref;
  std::cout << "Prop: Norm of difference GP vs Xconj for unit source special case (expect 0): " << norm2(diff) << std::endl;
  assert(norm2(diff) < 1e-8);
#endif


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
