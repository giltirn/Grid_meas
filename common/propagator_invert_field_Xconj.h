#pragma once

#include "propagator_invert.h"

namespace GridMeas{
  using namespace Grid; 

  //Invert the propagator against a field using the X-conjugate action under the hood
  //The computational cost is the same as the regular G-parity action but it may be faster due to a more optimal implementation of the 1f actions
  //and the use of 1f eigenvectors internally
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertFieldXconj(FermionFieldD &sol, const FermionFieldD &src, FermionActionD &xconj_action_d, FermionActionF &xconj_action_f, 
			  const MixedCGargs &args,
			  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    conformable(src.Grid(), sol.Grid());
    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermionField1fD> hermop_d(xconj_action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionField1fF> hermop_f(xconj_action_f);

    LinearFunction<FermionField1fD>* mcg = MixedCGfactory(args, xconj_action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
    LinearFunctionWrapper<FermionField1fD> mcg_wrap(*mcg);

    LinearFunction<FermionField1fD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermionField1fD>(*evals,*evecs);

    //ConjugateGradient<FermionField> CG(tol,10000);
    SchurRedBlackDiagMooeeSolve<FermionField1fD> solver(mcg_wrap);
  
    GridBase* FGridD = xconj_action_d.FermionGrid();
 
    //Extract the X-conjugate component vectors
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;


    FermionField1fD phi1 = PeekIndex<GparityFlavourIndex>(src,0);
    FermionField1fD phi2 = PeekIndex<GparityFlavourIndex>(src,1);
    FermionField1fD rho = 0.5*( phi1 + X*conjugate(phi2));
    FermionField1fD tau = ComplexD(0,0.5)*( phi1 - X*conjugate(phi2));
    
    FermionField1fD* srcs[2] = {&rho,&tau};
    sol = Zero();

    FermionField1fD src_5d(FGridD);
    FermionField1fD sol_5d(FGridD);    
    FermionField1fD sol_4d(src.Grid());
    FermionFieldD sol_4d_2f(src.Grid());
 
    for(int i=0;i<2;i++){
      xconj_action_d.ImportPhysicalFermionSource(*srcs[i], src_5d);

      guesser != nullptr ? 
	solver(xconj_action_d, src_5d, sol_5d, *guesser) 
	: 
	solver(xconj_action_d, src_5d, sol_5d);
      
      xconj_action_d.ExportPhysicalFermionSolution(sol_5d, sol_4d);
      
      PokeIndex<GparityFlavourIndex>(sol_4d_2f, sol_4d, 0);
      sol_4d = -(X*conjugate(sol_4d));
      PokeIndex<GparityFlavourIndex>(sol_4d_2f, sol_4d, 1);

      if(i==0) sol = sol_4d_2f;
      else sol = sol + ComplexD(0,-1)*sol_4d_2f;
    }    

    if(guesser) delete guesser;
    delete mcg;
  }

  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertFieldXconj(FermionFieldD &sol, const FermionFieldD &src, FermionActionD &xconj_action_d, FermionActionF &xconj_action_f, 
				 const MixedCGargs &args){
    mixedPrecInvertFieldXconj(sol,src,xconj_action_d,xconj_action_f,args,(std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }
};
