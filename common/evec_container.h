#pragma once

#include "propagator_invert.h"
#include "propagator_invert_Xconj.h"
#include "utils.h"
#include "lanczos.h"

//A container class for generating and utilizing eigenvectors that works for both X-conjugate and G-parity actions
namespace GridMeas{
  using namespace Grid;

  struct EvecContainerOpts{
    bool load_evecs;
    std::string load_evecs_stub;
    std::string load_evals_stub;

    bool save_evecs;
    std::string save_evecs_stub;
    std::string save_evals_stub;

    EvecContainerOpts(): load_evecs(false), save_evecs(false){}
  };

  template<typename _LanczosAction>
  struct EvecContainer{
    typedef _LanczosAction LanczosAction;
    typedef typename LanczosAction::FermionField EvecFieldType;
  
    std::vector<RealD> eval;
    std::vector<EvecFieldType> evec;

    EvecContainer(){}

    void clear(){
      eval.clear();
      std::vector<EvecFieldType>().swap(evec);
    }

    void generate(const LanczosParameters &lanc_arg, int traj, 
		  LanczosAction &action, 
		  const typename LanczosAction::GaugeField &U, GridParallelRNG &pRNG, 
		  const EvecContainerOpts &opts= EvecContainerOpts() ){
      GridCartesian* FGrid = (GridCartesian*)action.FermionGrid();
      GridRedBlackCartesian* FrbGrid = (GridRedBlackCartesian*)action.FermionRedBlackGrid();

      clear();
      if(opts.load_evecs){
	readEigenvalues(eval, evec, FrbGrid, opts.load_evals_stub, opts.load_evecs_stub, traj);
      }else{
	computeEigenvalues<LanczosAction, EvecFieldType>(eval, evec, lanc_arg, FGrid, FrbGrid, U, action, pRNG);
      }
      if(opts.save_evecs) saveEigenvalues(eval, evec, opts.save_evals_stub, opts.save_evecs_stub, traj);
    }

    //5D solutions
    template<typename ActionD, typename ActionF>
    void mixedPrecInvert5D(typename ActionD::PropagatorField &prop5d, 
				    const typename ActionD::PropagatorField &msrc, ActionD &action_d, ActionF &action_f, 
				    const MixedCGargs &args){
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      prop5d = GridMeas::mixedPrecInvert5D(msrc, action_d, action_f, args, eval_ptr, evec_ptr);
    }
    template<typename ActionD, typename ActionF>
    void splitGridMixedPrecInvert5D(typename ActionD::PropagatorField &msol5d,
				    const typename ActionD::PropagatorField &msrc,
				    ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
				    const MixedCGargs &args){
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::splitGridMixedPrecInvert5D(msol5d, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval_ptr, evec_ptr);
    }

    //4D solutions
    template<typename ActionD, typename ActionF>
    void mixedPrecInvertWithMidProp(typename ActionD::PropagatorField &prop, typename ActionD::PropagatorField &midprop, 
				    const typename ActionD::PropagatorField &msrc, ActionD &action_d, ActionF &action_f, 
				    const MixedCGargs &args){
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
    }  

    template<typename ActionD, typename ActionF>
    void splitGridMixedPrecInvertWithMidProp(typename ActionD::PropagatorField &msol, typename ActionD::PropagatorField &msol_mid,
					     const typename ActionD::PropagatorField &msrc,
					     ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
					     const MixedCGargs &args){
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval_ptr, evec_ptr);
    }

    //These versions use the X-conjugate actions to compute G-parity evecs (2f)
    template<typename ActionD, typename ActionF, typename std::enable_if<isXconjAction<ActionD>::value, int>::type = 0  >
    void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				    const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
				    const MixedCGargs &args){
      std::cout << GridLogMessage <<"Inverting with X-conjugate action" << std::endl;
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::mixedPrecInvertWithMidPropXconj(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
    }
    template<typename ActionD, typename ActionF, typename std::enable_if<isXconjAction<ActionD>::value, int>::type = 0 >
    void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					     const LatticeSCFmatrixD &msrc,
					     ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
						    const MixedCGargs &args){
      std::cout << GridLogMessage << "Inverting with X-conjugate action" << std::endl;
      assert(0);
    }
  };





}
