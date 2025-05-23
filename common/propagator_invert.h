#pragma once

#include "defines.h"
#include "grids.h"
#include "guesser.h"
#include "field_utils.h"

namespace GridMeas{
  using namespace Grid;

  template<class FieldD>
  class LinearFunctionWrapper: public OperatorFunction<FieldD> {
    LinearFunction<FieldD> &mcg;

  public:
    LinearFunctionWrapper(LinearFunction<FieldD> &mcg): mcg(mcg){}

    void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &in, FieldD &out){
      mcg(in, out);
    }
  }; 

  GRID_SERIALIZABLE_ENUM(MixedCGalgorithm, undef, RestartedCG,1, ReliableUpdateCG,2);
  
  struct MixedCGargs: Serializable{
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedCGargs,
				    MixedCGalgorithm, algorithm,
				    double, tolerance,
				    double, restartedcg_inner_tol,
				    double, relupcg_delta);
    MixedCGargs(){
      algorithm = MixedCGalgorithm::RestartedCG;
      tolerance = 1e-8;
      restartedcg_inner_tol = 1e-5;
      relupcg_delta = 0.5;
    }
  };
      
  template<typename FieldD, typename FieldF>
  LinearFunction<FieldD>* MixedCGfactory(const MixedCGargs &args,
					 GridBase* sp_grid, 
					 LinearOperatorBase<FieldF> &Linop_f, 
					 LinearOperatorBase<FieldD> &Linop_d){
    if(args.algorithm == MixedCGalgorithm::RestartedCG){
      MixedPrecisionConjugateGradient<FieldD, FieldF>* out = new MixedPrecisionConjugateGradient<FieldD, FieldF>(args.tolerance, 10000,10000, sp_grid, Linop_f, Linop_d);
      out->InnerTolerance = args.restartedcg_inner_tol;
      return out;
    }else if(args.algorithm == MixedCGalgorithm::ReliableUpdateCG){
      return new ConjugateGradientReliableUpdate<FieldD, FieldF>(args.tolerance, 10000, args.relupcg_delta, sp_grid, Linop_f, Linop_d);
    }else{
      std::cout << "MixedCGfactory: unknown algorithm" << std::endl;
      assert(0);
    }
  }				 

  //General mixed prec inverter with optional midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertGen(typename FermionActionD::PropagatorField &prop, 
			  typename FermionActionD::PropagatorField &midprop, 
			  const typename FermionActionD::PropagatorField &msrc, 
			  FermionActionD &action_d, FermionActionF &action_f,
			  const MixedCGargs &args,
 			  bool do_midprop,
			  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typedef typename FermionActionD::FermionField FermTypeD;
    typedef typename FermionActionF::FermionField FermTypeF;

    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermTypeD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermTypeF> hermop_f(action_f);
  
    LinearFunction<FermTypeD>* mcg = MixedCGfactory(args, action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
    LinearFunctionWrapper<FermTypeD> mcg_wrap(*mcg);

    LinearFunction<FermTypeD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermTypeD>(*evals,*evecs);

    SchurRedBlackDiagMooeeSolve<FermTypeD> solver(mcg_wrap);
  
    GridBase* FGridD = action_d.FermionGrid();

    FermTypeD src_5d(FGridD);
    FermTypeD sol_5d(FGridD);
  
    FermTypeD sol_4d(msrc.Grid());
    conformable(msrc.Grid(), prop.Grid());
    if(do_midprop) conformable(msrc.Grid(), midprop.Grid());

    //Columns of the matrix are the source vectors
    typedef columnOps<FermTypeD> cop;
    int Ncols = cop::Ncols;

    for(int col=0;col<Ncols;col++){
      std::cout << GridLogMessage << "Starting column inversion " << (col+1) << "/" << Ncols << std::endl;
      std::cout << GridLogMessage << "Extracting source matrix column" << std::endl;
      FermTypeD src_4d = cop::extractColumn(msrc, col);
      std::cout << GridLogMessage << "Generating 5D source" << std::endl;
      action_d.ImportPhysicalFermionSource(src_4d, src_5d);
      
      std::cout << GridLogMessage << "Inverting" << std::endl;
      guesser != nullptr ? 
	solver(action_d, src_5d, sol_5d, *guesser) 
	: 
	solver(action_d, src_5d, sol_5d);
      
      std::cout << GridLogMessage << "Generating 4D solution" << std::endl;
      action_d.ExportPhysicalFermionSolution(sol_5d, sol_4d);
      
      std::cout << GridLogMessage << "Extracting solution matrix column" << std::endl;
      cop::insertColumn(prop, sol_4d, col);

      if(do_midprop){
	std::cout << GridLogMessage << "Generating 4D midpoint solution" << std::endl;
	sol_4d = extractMidProp(sol_5d, action_d);
	
	std::cout << GridLogMessage << "Extracting solution midpoint matrix column" << std::endl;
	cop::insertColumn(midprop, sol_4d, col);	  
      }
    }

    if(guesser) delete guesser;
    delete mcg;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
  }

  //No midprop, with evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  typename FermionActionD::PropagatorField mixedPrecInvert(const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
							   const MixedCGargs &args,
							   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typename FermionActionD::PropagatorField tmp(msrc.Grid()), prop(msrc.Grid());
    mixedPrecInvertGen(prop,tmp,msrc,action_d,action_f,args,false,evals,evecs);
    return prop;
  }
  //No midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  typename FermionActionD::PropagatorField mixedPrecInvert(const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
							   const MixedCGargs &args){
    return mixedPrecInvert(msrc,action_d,action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //With midprop and evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &prop, typename FermionActionD::PropagatorField &midprop, 
				  const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				  const MixedCGargs &args,
				  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    mixedPrecInvertGen(prop,midprop,msrc,action_d,action_f,args,true,evals,evecs);
  }    
  //With midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &prop, typename FermionActionD::PropagatorField &midprop, 
				  const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				  const MixedCGargs &args){
    mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, args, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }

  //For Nsrc sources, the Nsrc * ( 24(G-parity) or 12(periodic) ) inversions will be farmed out to the subgrids using the split-grid method
  //subgrid_action_f must be defined on the subgrids
  //subgrids_f must be the single precision subgrids
  //if do_midprop = true the msol_midprop output will be populated with the midpoint solutions
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertGen(std::vector<typename FermionActionD::PropagatorField> &msol, std::vector<typename FermionActionD::PropagatorField> &msol_midprop,
				   const std::vector<typename FermionActionD::PropagatorField> &msrc,
				   FermionActionD &action_d,
				   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				   const MixedCGargs &args,
				   bool do_midprop,
				   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typedef typename FermionActionD::FermionField FermTypeD;
    typedef typename FermionActionF::FermionField FermTypeF;
    int Nsrc = msrc.size();
    
    //Setup outputs
    msol.resize(Nsrc, msrc[0].Grid());
    msol_midprop.resize( do_midprop ? Nsrc : 0, msrc[0].Grid() );
    
    //Compute the number of subgrids
    GridCartesian* Ugrid_full = (GridCartesian*)action_d.GaugeGrid();
    GridCartesian* Ugrid_sub = (GridCartesian*)subgrid_action_d.GaugeGrid();
    int Nsubgrids = 1;
    for(int i=0;i<Nd;i++)
      Nsubgrids *= Ugrid_full->_processors[i]/Ugrid_sub->_processors[i];

    std::cout << GridLogMessage << "Starting split Grid inversion with " << Nsubgrids << " subgrids" << std::endl;
    
    //Setup deflation
    LinearFunction<FermTypeD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermTypeD>(*evals,*evecs);

    bool deflate = guesser != nullptr;
    
    //Setup the subgrid solver
    SchurDiagMooeeOperator<FermionActionD,FermTypeD> hermop_d(subgrid_action_d);
    SchurDiagMooeeOperator<FermionActionF,FermTypeF> hermop_f(subgrid_action_f);

    LinearFunction<FermTypeD>* mcg = MixedCGfactory(args, subgrid_action_f.FermionRedBlackGrid(),  hermop_f, hermop_d);
    LinearFunctionWrapper<FermTypeD> mcg_wrap(*mcg);  
    SchurRedBlackDiagMooeeSolve<FermTypeD> solver(mcg_wrap);
    
    //Setup scratch space for the subgrid solves
    std::vector<FermTypeD> split_sources_5d_full(Nsubgrids, action_d.FermionGrid());
    std::vector<FermTypeD> split_solutions_5d_full(Nsubgrids, action_d.FermionGrid());
    FermTypeD split_source_5d(subgrid_action_d.FermionGrid()); //the subgrid source
    FermTypeD split_solution_5d(subgrid_action_d.FermionGrid());

    //Solve in block    
    typedef columnOps<FermTypeD> cop;
    int Ncols = cop::Ncols;
    int Ninvert = Ncols*Nsrc;
    int Nblocks = (Ninvert + Nsubgrids -1)/Nsubgrids;
    std::cout << GridLogMessage << "Doing split Grid inversion of " << Ninvert << " source vectors in " << Nblocks << " blocks" << std::endl;
    for(int b=0;b<Nblocks;b++){
      int istart = b*Nsubgrids;
      int ilessthan = std::min(istart + Nsubgrids, Ninvert);
      int sidx=0; //split source/soln index
      std::cout << GridLogMessage << "Generating 5d sources for block " << b << std::endl;
      for(int i=istart;i<ilessthan;i++){ //i = c + Nc *( s + Ns*( f + Ngp * src ) ) 
	int pcol = i % Ncols;
	int src = i / Ncols;

	FermTypeD src_4d = cop::extractColumn(msrc[src], pcol);
	action_d.ImportPhysicalFermionSource(src_4d, split_sources_5d_full[sidx++]);
      }	
      while(sidx < Nsubgrids)
	split_sources_5d_full[sidx++] = split_sources_5d_full[0]; //give the wasted subgrids something to work on

      std::cout << GridLogMessage << "Splitting 5d sources" << std::endl;
      Grid_split(split_sources_5d_full, split_source_5d);

      //Setup guess
      if(deflate){
	std::cout << GridLogMessage << "Deflating 5d sources" << std::endl;
	sidx=0;
	for(int i=istart;i<ilessthan;i++){
	  (*guesser)(split_sources_5d_full[sidx], split_solutions_5d_full[sidx]);
	  ++sidx;
	}
	while(sidx < Nsubgrids)
	  split_solutions_5d_full[sidx++] = split_solutions_5d_full[0];

	std::cout << GridLogMessage << "Splitting guesses" << std::endl;
	Grid_split(split_solutions_5d_full, split_solution_5d);
      }else{
	split_solution_5d = Zero();
      }
      
      //Subgrid solve
      std::cout << GridLogMessage << "Doing subgrid solve for block " << b << std::endl;
      solver(subgrid_action_d, split_source_5d, split_solution_5d);

      std::cout << GridLogMessage << "Unsplitting solutions" << std::endl;
      Grid_unsplit(split_solutions_5d_full, split_solution_5d);

      std::cout << GridLogMessage << "Extracting 4d solutions" << std::endl;
      FermTypeD sol_4d(msrc[0].Grid());
      sidx=0;
      for(int i=istart;i<ilessthan;i++){
	int pcol = i % Ncols;
	int src = i / Ncols;
	
	action_d.ExportPhysicalFermionSolution(split_solutions_5d_full[sidx], sol_4d);
	cop::insertColumn(msol[src], sol_4d, pcol);
	
	if(do_midprop){
	  sol_4d = extractMidProp(split_solutions_5d_full[sidx], action_d);
	  cop::insertColumn(msol_midprop[src], sol_4d, pcol);	  
	}
	++sidx;
      }
    }    
    if(deflate) delete guesser;	
    delete mcg;
  }

  //Multi-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert(std::vector<typename FermionActionD::PropagatorField> &msol,
				const std::vector<typename FermionActionD::PropagatorField> &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				const MixedCGargs &args,
				std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<typename FermionActionD::PropagatorField> tmp;
    splitGridMixedPrecInvertGen(msol, tmp, msrc, action_d, subgrid_action_d, subgrid_action_f, args, false, evals, evecs);
  }
  //Multi-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvert(std::vector<typename FermionActionD::PropagatorField> &msol,
				const std::vector<typename FermionActionD::PropagatorField> &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				const MixedCGargs &args){
    splitGridMixedPrecInvert(msol, msrc, action_d, subgrid_action_d, subgrid_action_f, args, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //Single-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert(typename FermionActionD::PropagatorField &msol,
				const typename FermionActionD::PropagatorField &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				const MixedCGargs &args,
				std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    
    std::vector<typename FermionActionD::PropagatorField> tmp_msrc(1, msrc);
    std::vector<typename FermionActionD::PropagatorField> tmp_msol;
    std::vector<typename FermionActionD::PropagatorField> tmp_msol_mid;
    splitGridMixedPrecInvertGen(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, args, false, evals, evecs);
    msol = tmp_msol[0];
  }
  //Single-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvert(typename FermionActionD::PropagatorField &msol,
				const typename FermionActionD::PropagatorField &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				const MixedCGargs &args){
    splitGridMixedPrecInvert(msol,msrc,action_d,subgrid_action_d,subgrid_action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }
				


  //Multi-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidProp(std::vector<typename FermionActionD::PropagatorField> &msol, std::vector<typename FermionActionD::PropagatorField> &msol_mid,
					   const std::vector<typename FermionActionD::PropagatorField> &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   const MixedCGargs &args,
					   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    splitGridMixedPrecInvertGen(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, true, evals, evecs);
  }
  //Multi-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidProp(std::vector<typename FermionActionD::PropagatorField> &msol, std::vector<typename FermionActionD::PropagatorField> &msol_mid,
					   const std::vector<typename FermionActionD::PropagatorField> &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   const MixedCGargs &args){
    splitGridMixedPrecInvertWithMidProp(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //Single-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &msol, typename FermionActionD::PropagatorField &msol_mid,
					   const typename FermionActionD::PropagatorField &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   const MixedCGargs &args,
					   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<typename FermionActionD::PropagatorField> tmp_msrc(1, msrc);
    std::vector<typename FermionActionD::PropagatorField> tmp_msol;
    std::vector<typename FermionActionD::PropagatorField> tmp_msol_mid;        
    splitGridMixedPrecInvertGen(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, args, true, evals, evecs);
    msol = tmp_msol[0];
    msol_mid = tmp_msol_mid[0];
  }
  //Single-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &msol, typename FermionActionD::PropagatorField &msol_mid,
					   const typename FermionActionD::PropagatorField &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   const MixedCGargs &args){
    splitGridMixedPrecInvertWithMidProp(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,args, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }




  

};
