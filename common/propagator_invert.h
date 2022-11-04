#pragma once

#include "defines.h"
#include "grids.h"
#include "guesser.h"
#include "field_utils.h"

namespace GridMeas{
  using namespace Grid;

  template<class FieldD, class FieldF>
  class MixedCGwrapper: public OperatorFunction<FieldD> {
    MixedPrecisionConjugateGradient<FieldD, FieldF> &mcg;

  public:
    MixedCGwrapper(MixedPrecisionConjugateGradient<FieldD, FieldF> &mcg): mcg(mcg){}

    void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &in, FieldD &out){
      mcg(in, out);
    }
  }; 

  //General mixed prec inverter with optional midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertGen(typename FermionActionD::PropagatorField &prop, 
			  typename FermionActionD::PropagatorField &midprop, 
			  const typename FermionActionD::PropagatorField &msrc, 
			  FermionActionD &action_d, FermionActionF &action_f, 
			  double tol, double inner_tol,
			  bool do_midprop,
			  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typedef typename FermionActionD::FermionField FermTypeD;
    typedef typename FermionActionF::FermionField FermTypeF;

    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermTypeD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermTypeF> hermop_f(action_f);
  
    MixedPrecisionConjugateGradient<FermTypeD, FermTypeF> mcg(tol, 10000,10000, action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper<FermTypeD, FermTypeF> mcg_wrap(mcg);

    LinearFunction<FermTypeD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermTypeD>(*evals,*evecs);

    //ConjugateGradient<FermionField> CG(tol,10000);
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
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
  }

  //No midprop, with evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  typename FermionActionD::PropagatorField mixedPrecInvert(typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
							   double tol, double inner_tol,
							   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typename FermionActionD::PropagatorField tmp(msrc.Grid()), prop(msrc.Grid());
    mixedPrecInvertGen(prop,tmp,msrc,action_d,action_f,tol,inner_tol,false,evals,evecs);
    return prop;
  }
  //No midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  typename FermionActionD::PropagatorField mixedPrecInvert(const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
							   double tol, double inner_tol){
    return mixedPrecInvert(msrc,action_d,action_f,tol,inner_tol,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //With midprop and evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &prop, typename FermionActionD::PropagatorField &midprop, 
				  const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				  double tol, double inner_tol,
				  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    mixedPrecInvertGen(prop,midprop,msrc,action_d,action_f,tol,inner_tol,true,evals,evecs);
  }    
  //With midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertWithMidProp(typename FermionActionD::PropagatorField &prop, typename FermionActionD::PropagatorField &midprop, 
				  const typename FermionActionD::PropagatorField &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				  double tol, double inner_tol){
    mixedPrecInvertWithMidProp(prop, midprop, msrc, action_d, action_f, tol, inner_tol, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }

  //For Nsrc sources, the 24*Nsrc inversions will be farmed out to the subgrids using the split-grid method
  //subgrid_action_f must be defined on the subgrids
  //subgrids_f must be the single precision subgrids
  //if do_midprop = true the msol_midprop output will be populated with the midpoint solutions
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertGen(std::vector<LatticeSCFmatrixD> &msol, std::vector<LatticeSCFmatrixD> &msol_midprop,
				   const std::vector<LatticeSCFmatrixD> &msrc,
				   FermionActionD &action_d,
				   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				   double tol, double inner_tol,
				   bool do_midprop,
				   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
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
    LinearFunction<FermionFieldD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermionFieldD>(*evals,*evecs);

    bool deflate = guesser != nullptr;
    
    //Setup the subgrid solver
    SchurDiagMooeeOperator<FermionActionD,FermionFieldD> hermop_d(subgrid_action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionFieldF> hermop_f(subgrid_action_f);

    MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> mcg(tol, 10000,10000, subgrid_action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper<FermionFieldD, FermionFieldF> mcg_wrap(mcg);  
    SchurRedBlackDiagMooeeSolve<FermionFieldD> solver(mcg_wrap);
    
    //Setup scratch space for the subgrid solves
    std::vector<FermionFieldD> split_sources_5d_full(Nsubgrids, action_d.FermionGrid());
    std::vector<FermionFieldD> split_solutions_5d_full(Nsubgrids, action_d.FermionGrid());
    FermionFieldD split_source_5d(subgrid_action_d.FermionGrid()); //the subgrid source
    FermionFieldD split_solution_5d(subgrid_action_d.FermionGrid());

    //Solve in blocks
    int Ninvert = Ngp*Ns*Nc*Nsrc;
    int Nblocks = (Ninvert + Nsubgrids -1)/Nsubgrids;
    std::cout << GridLogMessage << "Doing split Grid inversion of " << Ninvert << " source vectors in " << Nblocks << " blocks" << std::endl;
    for(int b=0;b<Nblocks;b++){
      int istart = b*Nsubgrids;
      int ilessthan = std::min(istart + Nsubgrids, Ninvert);
      int sidx=0; //split source/soln index
      std::cout << GridLogMessage << "Generating 5d sources for block " << b << std::endl;
      for(int i=istart;i<ilessthan;i++){ //i = c + Nc *( s + Ns*( f + Ngp * src ) ) 
	int rem = i;
	int c = rem % Nc; rem /= Nc;
	int s = rem % Ns; rem /= Ns;
	int f = rem % Ngp; rem /= Ngp;
	int src = rem;

	FermionFieldD src_4d = extractColumn(msrc[src], f,s,c);
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
      FermionFieldD sol_4d(msrc[0].Grid());
      sidx=0;
      for(int i=istart;i<ilessthan;i++){
	int rem = i;
	int c = rem % Nc; rem /= Nc;
	int s = rem % Ns; rem /= Ns;
	int f = rem % Ngp; rem /= Ngp;
	int src = rem;
	
	action_d.ExportPhysicalFermionSolution(split_solutions_5d_full[sidx], sol_4d);
	insertColumn(msol[src], sol_4d, f,s,c);
	
	if(do_midprop){
	  sol_4d = extractMidProp(split_solutions_5d_full[sidx], action_d);
	  insertColumn(msol_midprop[src], sol_4d, f,s,c);	  
	}
	++sidx;
      }
    }
    if(deflate) delete guesser;	
  }

  //Multi-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert(std::vector<LatticeSCFmatrixD> &msol,
				const std::vector<LatticeSCFmatrixD> &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				double tol, double inner_tol,
				std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<LatticeSCFmatrixD> tmp;
    splitGridMixedPrecInvertGen(msol, tmp, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, false, evals, evecs);
  }
  //Multi-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvert(std::vector<LatticeSCFmatrixD> &msol,
				const std::vector<LatticeSCFmatrixD> &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				double tol, double inner_tol){
    splitGridMixedPrecInvert(msol, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
  }


  //Single-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert(LatticeSCFmatrixD &msol,
				const LatticeSCFmatrixD &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				double tol, double inner_tol,
				std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    
    std::vector<LatticeSCFmatrixD> tmp_msrc(1, msrc);
    std::vector<LatticeSCFmatrixD> tmp_msol;
    std::vector<LatticeSCFmatrixD> tmp_msol_mid;
    splitGridMixedPrecInvertGen(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, false, evals, evecs);
    msol = tmp_msol[0];
  }
  //Single-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvert(LatticeSCFmatrixD &msol,
				const LatticeSCFmatrixD &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				double tol, double inner_tol){
    splitGridMixedPrecInvert(msol,msrc,action_d,subgrid_action_d,subgrid_action_f,tol,inner_tol,(std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
  }
				


  //Multi-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidProp(std::vector<LatticeSCFmatrixD> &msol, std::vector<LatticeSCFmatrixD> &msol_mid,
					   const std::vector<LatticeSCFmatrixD> &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   double tol, double inner_tol,
					   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    splitGridMixedPrecInvertGen(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, true, evals, evecs);
  }
  //Multi-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidProp(std::vector<LatticeSCFmatrixD> &msol, std::vector<LatticeSCFmatrixD> &msol_mid,
					   const std::vector<LatticeSCFmatrixD> &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   double tol, double inner_tol){
    splitGridMixedPrecInvertWithMidProp(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,tol,inner_tol,(std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
  }


  //Single-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					   const LatticeSCFmatrixD &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   double tol, double inner_tol,
					   std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<LatticeSCFmatrixD> tmp_msrc(1, msrc);
    std::vector<LatticeSCFmatrixD> tmp_msol;
    std::vector<LatticeSCFmatrixD> tmp_msol_mid;        
    splitGridMixedPrecInvertGen(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, tol, inner_tol, true, evals, evecs);
    msol = tmp_msol[0];
    msol_mid = tmp_msol_mid[0];
  }
  //Single-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidProp(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
					   const LatticeSCFmatrixD &msrc,
					   FermionActionD &action_d,
					   FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
					   double tol, double inner_tol){
    splitGridMixedPrecInvertWithMidProp(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,tol,inner_tol, (std::vector<Real> const*)nullptr, (std::vector<FermionFieldD> const *)nullptr);
  }




  

};
