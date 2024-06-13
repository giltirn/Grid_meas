#pragma once

#include "propagator_invert.h"

namespace GridMeas{
  using namespace Grid; 


  //Invert the propagator against a field using the X-conjugate action under the hood
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertField(typename FermionActionD::FermionField &sol, const typename FermionActionD::FermionField &src, FermionActionD &action_d, FermionActionF &action_f, 
			  const MixedCGargs &args,
			  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typedef typename FermionActionD::FermionField FermionFieldD;
    typedef typename FermionActionF::FermionField FermionFieldF;

    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermionFieldD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionFieldF> hermop_f(action_f);
  
    LinearFunction<FermionFieldD>* mcg = MixedCGfactory(args, action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
    LinearFunctionWrapper<FermionFieldD> mcg_wrap(*mcg);

    LinearFunction<FermionFieldD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermionFieldD>(*evals,*evecs);

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
    delete mcg;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
  }

  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertField(typename FermionActionD::FermionField &sol, const typename FermionActionD::FermionField &src, FermionActionD &action_d, FermionActionF &action_f, 
			    const MixedCGargs &args){
    mixedPrecInvertField(sol,src,action_d,action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }

  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvert5dField(typename FermionActionD::FermionField &sol, const typename FermionActionD::FermionField &src, FermionActionD &action_d, FermionActionF &action_f, 
			      const MixedCGargs &args){
    typedef typename FermionActionD::FermionField FermionFieldD;
    typedef typename FermionActionF::FermionField FermionFieldF;

    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermionFieldD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionFieldF> hermop_f(action_f);
  
    LinearFunction<FermionFieldD>* mcg = MixedCGfactory(args, action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
    
    (*mcg)(src,sol);
  }

    
};
