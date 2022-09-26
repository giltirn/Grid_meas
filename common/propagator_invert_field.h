#pragma once

#include "propagator_invert.h"

namespace GridMeas{
  using namespace Grid; 

 //General mixed prec inverter with optional midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertField(FermionFieldD &sol, const FermionFieldD &src, FermionActionD &action_d, FermionActionF &action_f, 
			  double tol, double inner_tol,
			  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermionFieldD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionFieldF> hermop_f(action_f);
  
    MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> mcg(tol, 10000,10000, action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper<FermionFieldD, FermionFieldF> mcg_wrap(mcg);

    LinearFunction<FermionFieldD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermionFieldD>(*evals,*evecs);

    //ConjugateGradient<FermionField> CG(tol,10000);
    SchurRedBlackDiagMooeeSolve<FermionFieldD> solver(mcg_wrap);
  
    GridBase* FGridD = action_d.FermionGrid();

    FermionFieldD src_5d(FGridD);
    FermionFieldD sol_5d(FGridD);
  
    FermionFieldD sol_4d(src.Grid());
    conformable(src.Grid(), sol.Grid());
    
    action_d.ImportPhysicalFermionSource(src, src_5d);
    
    std::cout << GridLogMessage << "Inverting" << std::endl;
    guesser != nullptr ? 
      solver(action_d, src_5d, sol_5d, *guesser) 
      : 
      solver(action_d, src_5d, sol_5d);
	
    std::cout << GridLogMessage << "Generating 4D solution" << std::endl;
    action_d.ExportPhysicalFermionSolution(sol_5d, sol);
    
    if(guesser) delete guesser;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
  }


};
