#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //From the matrices extract the column with the given flavor, spin and color index
  FermionFieldD extractColumn(const LatticeSCFmatrixD &from, const int fs, const int ss, const int cs){
    FermionFieldD into(from.Grid());
    {
      size_t Nsimd = FermionFieldD::vector_type::Nsimd();
      autoView( from_v, from, AcceleratorRead);
      autoView( into_v, into, AcceleratorWrite);
      accelerator_for( i, from_v.size(), Nsimd, {
	  auto site_matrix = from_v(i);
	  auto site_into = into_v(i);
	  for(int fr=0;fr<Ngp;fr++)
	    for(int sr=0;sr<Ns;sr++)
	      for(int cr=0;cr<Nc;cr++)
		site_into(fr)(sr)(cr) = site_matrix(fr,fs)(sr,ss)(cr,cs);
	  coalescedWrite(into_v[i], site_into);
	});
    }
    return into;
  }

  void insertColumn(LatticeSCFmatrixD &into, const FermionFieldD &from, const int fs, const int ss, const int cs){
    size_t Nsimd = FermionFieldD::vector_type::Nsimd();
    autoView( from_v, from, AcceleratorRead);
    autoView( into_v, into, AcceleratorWrite);
    accelerator_for( i, from_v.size(), Nsimd, {
	auto site_spinor = from_v(i);
	auto site_into = into_v(i);
	for(int fr=0;fr<Ngp;fr++)
	  for(int sr=0;sr<Ns;sr++)
	    for(int cr=0;cr<Nc;cr++)
	      site_into(fr,fs)(sr,ss)(cr,cs) = site_spinor(fr)(sr)(cr);
	coalescedWrite(into_v[i], site_into);
      });
  }

  class MixedCGwrapper: public OperatorFunction<FermionFieldD> {
    MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> &mcg;

  public:
    MixedCGwrapper(MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> &mcg): mcg(mcg){}

    void operator() (LinearOperatorBase<FermionFieldD> &Linop, const FermionFieldD &in, FermionFieldD &out){
      mcg(in, out);
    }
  };
    
  template<typename FermionActionD, typename FermionActionF>
  LatticeSCFmatrixD mixedPrecInvert(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, double tol, double inner_tol,
				    std::vector<Real> const* evals = nullptr, std::vector<FermionFieldD> const * evecs = nullptr){
    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<FermionActionD,FermionFieldD> hermop_d(action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionFieldF> hermop_f(action_f);
  
    MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> mcg(tol, 10000,10000, action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper mcg_wrap(mcg);
  
    DeflatedGuesser<FermionFieldD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = new DeflatedGuesser<FermionFieldD>(*evecs, *evals);

    //ConjugateGradient<FermionField> CG(tol,10000);
    SchurRedBlackDiagMooeeSolve<FermionFieldD> solver(mcg_wrap);
  
    GridBase* FGridD = action_d.FermionGrid();

    FermionFieldD src_5d(FGridD);
    FermionFieldD sol_5d(FGridD);
  
    FermionFieldD sol_4d(msrc.Grid());
    LatticeSCFmatrixD msol(msrc.Grid());

    //Columns of the matrix are the source vectors
    for(int f=0;f<Ngp;f++){
      for(int s=0;s<Ns;s++){
	for(int c=0;c<Nc;c++){
	  std::cout << GridLogMessage << "Starting f="<< f << " s=" << s << " c=" << c << " inversion" << std::endl;
	  std::cout << GridLogMessage << "Extracting source matrix column" << std::endl;
	  FermionFieldD src_4d = extractColumn(msrc, f,s,c);
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
	  insertColumn(msol, sol_4d, f,s,c);
	}
      }
    }

    if(guesser) delete guesser;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;

    return msol;
  }


  //The midpoint propagator
  //\phi(x) = + P_L \psi(x,Ls/2) + P_R \psi(x,Ls/2-1) 
  template<typename Impl>
  typename Impl::FermionField extractMidProp(const typename Impl::FermionField &sol5d, CayleyFermion5D<Impl> &action_d){
    typedef typename Impl::FermionField Ferm;
    GridBase* UGrid = action_d.GaugeGrid();
    GridBase* FGrid = action_d.FermionGrid();
    
    Gamma G5(Gamma::Algebra::Gamma5);
    int Ls = action_d.Ls;
    Ferm p_plus (UGrid);			  
    Ferm p_minus(UGrid);
    Ferm p(UGrid);

    ExtractSlice(p_plus , sol5d, Ls/2-1 , 0);
    ExtractSlice(p_minus, sol5d, Ls/2   , 0);
    p_plus = p_plus + G5*p_plus;
    p_minus= p_minus - G5*p_minus;
    p=0.5*(p_plus+p_minus);
    return p;
  }


  void mixedPrecInvertWithMidProp(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				  const LatticeSCFmatrixD &msrc, CayleyFermion5D<GparityWilsonImplD> &action_d, CayleyFermion5D<GparityWilsonImplF> &action_f, 
				  double tol, double inner_tol,
				  std::vector<Real> const* evals = nullptr, std::vector<FermionFieldD> const * evecs = nullptr){
    std::cout << GridLogMessage << "Starting source inversion" << std::endl;
    SchurDiagMooeeOperator<CayleyFermion5D<GparityWilsonImplD>,FermionFieldD> hermop_d(action_d);
    SchurDiagMooeeOperator<CayleyFermion5D<GparityWilsonImplF>,FermionFieldF> hermop_f(action_f);
  
    MixedPrecisionConjugateGradient<FermionFieldD, FermionFieldF> mcg(tol, 10000,10000, action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper mcg_wrap(mcg);
  
    DeflatedGuesser<FermionFieldD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = new DeflatedGuesser<FermionFieldD>(*evecs, *evals);

    //ConjugateGradient<FermionField> CG(tol,10000);
    SchurRedBlackDiagMooeeSolve<FermionFieldD> solver(mcg_wrap);
  
    GridBase* FGridD = action_d.FermionGrid();

    FermionFieldD src_5d(FGridD);
    FermionFieldD sol_5d(FGridD);
  
    FermionFieldD sol_4d(msrc.Grid());
    conformable(msrc.Grid(), prop.Grid());
    conformable(msrc.Grid(), midprop.Grid());

    //Columns of the matrix are the source vectors
    for(int f=0;f<Ngp;f++){
      for(int s=0;s<Ns;s++){
	for(int c=0;c<Nc;c++){
	  std::cout << GridLogMessage << "Starting f="<< f << " s=" << s << " c=" << c << " inversion" << std::endl;
	  std::cout << GridLogMessage << "Extracting source matrix column" << std::endl;
	  FermionFieldD src_4d = extractColumn(msrc, f,s,c);
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
	  insertColumn(prop, sol_4d, f,s,c);

	  std::cout << GridLogMessage << "Generating 4D midpoint solution" << std::endl;
	  sol_4d = extractMidProp(sol_5d, action_d);

	  std::cout << GridLogMessage << "Extracting solution midpoint matrix column" << std::endl;
	  insertColumn(midprop, sol_4d, f,s,c);	  
	}
      }
    }

    if(guesser) delete guesser;
    std::cout << GridLogMessage << "Source inversion complete" << std::endl;
  }




};
