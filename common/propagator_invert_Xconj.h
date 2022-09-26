#pragma once

#include "propagator_invert.h"


namespace GridMeas{
  using namespace Grid;

  //0.5*(1+iX) in
  template<typename T>
  T mulPplusLeft(const T &in){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static ComplexD _I(0,1);
    T out = 0.5*(in + _I * (X*in));
    return out;
  }
  //0.5*(1+iX) in
  template<typename T>
  T mulPminusLeft(const T &in){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static ComplexD _mI(0,-1);
    T out = 0.5*(in + _mI * (X*in));
    return out;
  }

  //0.5*(1+iX) in
  template<typename T>
  T mulPplusRight(const T &in){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static ComplexD _I(0,1);
    T out = 0.5*(in + _I * (in*X));
    return out;
  }
  //0.5*(1+iX) in
  template<typename T>
  T mulPminusRight(const T &in){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static ComplexD _mI(0,-1);
    T out = 0.5*(in + _mI * (in*X));
    return out;
  }


  //Right multiply by U
  template<typename T>
    T mulURight(const T &in){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);

    T out = 0.5*(in + in*X) + 0.5*(in*sigma3 - (in*X)*sigma3);
    return out;
  }
  
  //Invert using the X-conjugate action
  //NOTE: the source *must* have flavor structure  
  // | phi   0   |
  // |  0   phi* |  
  //where phi is a spin-color matrix that commutes with X=C g5  
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertGenXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
			       const LatticeSCFmatrixD &msrc, FermionActionD &xconj_action_d, FermionActionF &xconj_action_f, 
			       double tol, double inner_tol,
			       bool do_midprop,
			       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);


    //Check the source has the right structure
    {
      auto msrc_00 = PeekIndex<GparityFlavourIndex>(msrc,0,0);
      auto msrc_11 = PeekIndex<GparityFlavourIndex>(msrc,1,1);
      decltype(msrc_00) diff = msrc_00 - conjugate(msrc_11);
      assert(norm2(diff) < 1e-12);
      auto msrc_01 = PeekIndex<GparityFlavourIndex>(msrc,0,1);
      auto msrc_10 = PeekIndex<GparityFlavourIndex>(msrc,1,0);   
      assert(norm2(msrc_01) < 1e-12);
      assert(norm2(msrc_10) < 1e-12);

      decltype(msrc_00) tmp = X*msrc_00 - msrc_00*X;
      assert(norm2(tmp) < 1e-12);
    }      

    GridBase* UGrid = xconj_action_d.GaugeGrid();
    GridBase* FGrid = xconj_action_d.FermionGrid();
    GridBase* FrbGrid = xconj_action_d.FermionRedBlackGrid();


    SchurDiagMooeeOperator<FermionActionD,FermionField1fD> hermop_d(xconj_action_d);
    SchurDiagMooeeOperator<FermionActionF,FermionField1fF> hermop_f(xconj_action_f);
  
    MixedPrecisionConjugateGradient<FermionField1fD, FermionField1fF> mcg(tol, 10000,10000, xconj_action_f.FermionRedBlackGrid(), hermop_f, hermop_d);
    mcg.InnerTolerance = inner_tol;
    MixedCGwrapper<FermionField1fD, FermionField1fF> mcg_wrap(mcg);
    SchurRedBlackDiagMooeeSolve<FermionField1fD> solver(mcg_wrap);

    // ConjugateGradient<FermionField1fD> cg(1e-08,10000);
    // SchurRedBlackDiagMooeeSolve<FermionField1fD> solver(cg);

    LinearFunction<FermionField1fD> *guesser = nullptr;
    if(evecs != nullptr && evals != nullptr)
      guesser = getGuesser<FermionField1fD>(*evals,*evecs);


    FermionField1fD tmp4d(UGrid);
    FermionField1fD src4d(UGrid);
    FermionField1fD src5d(FGrid), sol5d(FGrid);
    
    FermionFieldD tmp2f(UGrid);
    LatticeSCFmatrixD V_4d(UGrid), V_4d_mid(UGrid);
    
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
    	for(int pm=0;pm<2;pm++){
	  tmp2f = extractColumn(msrc, 0,s,c);
	  src4d = PeekIndex<GparityFlavourIndex>(tmp2f,0); //00 element
    	  src4d = pm == 0 ? mulPplusLeft(src4d) : mulPminusLeft(src4d);

    	  xconj_action_d.ImportPhysicalFermionSource(src4d, src5d);

	  guesser != nullptr ? 
	    solver(xconj_action_d, src5d, sol5d, *guesser) 
	    : 
	    solver(xconj_action_d, src5d, sol5d);

    	  xconj_action_d.ExportPhysicalFermionSolution(sol5d, tmp4d);

    	  //Generate 2f X-conjugate output
    	  PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 0);
    	  tmp4d = -(X*conjugate(tmp4d));
    	  PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 1);
    	  insertColumn(V_4d, tmp2f, pm,s,c);

	  if(do_midprop){
	    tmp4d = extractMidProp(sol5d, xconj_action_d);
	    PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 0);
	    tmp4d = -(X*conjugate(tmp4d));
	    PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 1);
	    insertColumn(V_4d_mid, tmp2f, pm,s,c);
	  }

    	}
      }
    }

    LatticeSCFmatrixD tmp = mulPplusRight(V_4d) + mulPminusRight(V_4d)*sigma1;  
    prop = mulURight(tmp);
    
    if(do_midprop){
      LatticeSCFmatrixD tmp = mulPplusRight(V_4d_mid) + mulPminusRight(V_4d_mid)*sigma1;  
      midprop = mulURight(tmp);
    }
    if(guesser) delete guesser;
  }
};
