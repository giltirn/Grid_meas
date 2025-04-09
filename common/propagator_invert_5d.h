#include "propagator_invert.h"

//Functionality for inverses with full 5D propagator output
namespace GridMeas{
  using namespace Grid;

  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  typename FermionActionD::PropagatorField mixedPrecInvert5D(const typename FermionActionD::PropagatorField &msrc, 
							     FermionActionD &action_d, FermionActionF &action_f,
							     const MixedCGargs &args,
							     std::vector<Real> const* evals = nullptr, std::vector<EvecFieldType> const * evecs = nullptr){
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
    typename FermionActionD::PropagatorField prop5d(FGridD);

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
           
      std::cout << GridLogMessage << "Extracting solution matrix column" << std::endl;
      cop::insertColumn(prop5d, sol_5d, col);
    }

    if(guesser) delete guesser;
    delete mcg;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
    return prop5d;
  }

  //For Nsrc sources, the Nsrc * ( 24(G-parity) or 12(periodic) ) inversions will be farmed out to the subgrids using the split-grid method
  //subgrid_action_f must be defined on the subgrids
  //subgrids_f must be the single precision subgrids
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert5D(std::vector<typename FermionActionD::PropagatorField> &msol,
				  const std::vector<typename FermionActionD::PropagatorField> &msrc,
				  FermionActionD &action_d,
				  FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				  const MixedCGargs &args,
				  std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    typedef typename FermionActionD::FermionField FermTypeD;
    typedef typename FermionActionF::FermionField FermTypeF;
    int Nsrc = msrc.size();
    
    assert(msrc[0].Grid() == action_d.GaugeGrid());

    //Setup outputs
    msol.resize(Nsrc, action_d.FermionGrid());
    
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
      sidx=0;
      for(int i=istart;i<ilessthan;i++){
	int pcol = i % Ncols;
	int src = i / Ncols;	
	cop::insertColumn(msol[src], split_solutions_5d_full[sidx], pcol);	
	++sidx;
      }
    }    
    if(deflate) delete guesser;	
    delete mcg;
  }
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvert5D(typename FermionActionD::PropagatorField &msol,
				const typename FermionActionD::PropagatorField &msrc,
				FermionActionD &action_d,
				FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				const MixedCGargs &args,
				std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    
    std::vector<typename FermionActionD::PropagatorField> tmp_msrc(1, msrc);
    std::vector<typename FermionActionD::PropagatorField> tmp_msol;
    splitGridMixedPrecInvert5D(tmp_msol, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, args, evals, evecs);
    msol = std::move(tmp_msol[0]);
  }
}
