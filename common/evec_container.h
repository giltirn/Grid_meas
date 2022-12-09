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

  template<bool isGparity>
  struct _action_call_wrapper{};

  template<>
  struct _action_call_wrapper<true>{
    template<typename ActionD, typename ActionF, typename EvecFieldType>
    static void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
					   const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
					   const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
      std::cout << GridLogMessage << "Inverting with G-parity action" << std::endl;
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
    }
    template<typename ActionD, typename ActionF, typename EvecFieldType>
    static void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						    const LatticeSCFmatrixD &msrc,
						    ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
						    const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
      std::cout << GridLogMessage << "Inverting with G-parity action" << std::endl;
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval_ptr, evec_ptr);
    }
  };
  template<>
  struct _action_call_wrapper<false>{
    template<typename ActionD, typename ActionF, typename EvecFieldType>
    static void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
					   const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
					   const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
      std::cout << GridLogMessage <<"Inverting with X-conjugate action" << std::endl;
      std::vector<RealD> const* eval_ptr = eval.size() ? &eval : (std::vector<RealD> const*)nullptr;
      std::vector<EvecFieldType> const* evec_ptr = eval.size() ? &evec : (std::vector<EvecFieldType> const*)nullptr;
      GridMeas::mixedPrecInvertWithMidPropXconj(prop, midprop, msrc, action_d, action_f, args, eval_ptr, evec_ptr);
    }
    template<typename ActionD, typename ActionF, typename EvecFieldType>
    static void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						    const LatticeSCFmatrixD &msrc,
						    ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
						    const MixedCGargs &args, const std::vector<RealD> &eval, const std::vector<EvecFieldType> &evec){
      std::cout << GridLogMessage << "Inverting with X-conjugate action" << std::endl;
      assert(0);
    }
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
  
    template<typename ActionD, typename ActionF>
    void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				    const LatticeSCFmatrixD &msrc, ActionD &action_d, ActionF &action_f, 
				    const MixedCGargs &args) const{
      _action_call_wrapper<ActionD::isGparity>::mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, eval, evec);
    }
    template<typename ActionD, typename ActionF>
    void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					     const LatticeSCFmatrixD &msrc,
					     ActionD &action_d, ActionD &subgrid_action_d, ActionF &subgrid_action_f,
					     const MixedCGargs &args) const{
      _action_call_wrapper<ActionD::isGparity>::splitGridMixedPrecInvertWithMidProp(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, eval, evec);
    }
  };





}
