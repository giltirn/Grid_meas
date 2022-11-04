#pragma once

#include "propagator_invert.h"
#include "utils.h"

namespace GridMeas{
  using namespace Grid;

  //0.5*(1+iX) in
  template<typename T>
  T mulPplusLeft(const T &in){
    static ComplexD _I(0,1);
    T out = 0.5*(in + _I * (Xmatrix()*in));
    return out;
  }
  //0.5*(1+iX) in
  template<typename T>
  T mulPminusLeft(const T &in){
    static ComplexD _mI(0,-1);
    T out = 0.5*(in + _mI * (Xmatrix()*in));
    return out;
  }

  //0.5*(1+iX) in
  template<typename T>
  T mulPplusRight(const T &in){
    static ComplexD _I(0,1);
    T out = 0.5*(in + _I * (in*Xmatrix()));
    return out;
  }
  //0.5*(1+iX) in
  template<typename T>
  T mulPminusRight(const T &in){
    static ComplexD _mI(0,-1);
    T out = 0.5*(in + _mI * (in*Xmatrix()));
    return out;
  }

  //H = |P+     P-   |
  //    |-XP-   -XP+ |  
  template<typename T>
  T mulHLeft(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    return mulPplusLeft(T(h1pi*in + h1mi*(sigma3*in))) + mulPminusLeft(T(h1mi*(sigma1*in) - h1mi*(sigma2*in)));
  }
  template<typename T>
  T mulHdagLeft(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    return mulPplusLeft(T(h1mi*in + h1pi*(sigma3*in))) + mulPminusLeft(T(h1pi*(sigma1*in) - h1pi*(sigma2*in)));
  }
  template<typename T>
  T mulHRight(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    return mulPplusRight(T(h1pi*in + h1mi*in*sigma3)) + mulPminusRight(T(h1mi*in*sigma1 - h1mi*in*sigma2));
  }
  template<typename T>
  inline T mulHdagRight(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    return mulPplusRight(T(h1mi*in + h1pi*in*sigma3)) + mulPminusRight(T(h1pi*in*sigma1 - h1pi*in*sigma2));
  }




  //Right multiply by U
  template<typename T>
    T mulURight(const T &in){
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);

    T out = 0.5*(in + in*Xmatrix()) + 0.5*(in*sigma3 - (in*Xmatrix())*sigma3);
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

      decltype(msrc_00) tmp = Xmatrix()*msrc_00 - msrc_00*Xmatrix();
      assert(norm2(tmp) < 1e-12);
    }      

    LatticeSCFmatrixD src_rotated = mulHRight(msrc); //the columns of this are X-conjugate
    LatticePropagatorD src_fcol(msrc.Grid());
    std::vector<LatticePropagatorD> prop_fcol(2, msrc.Grid()), midprop_fcol(2, msrc.Grid());

    for(int fcol=0;fcol<Ngp;fcol++){
      src_fcol = PeekIndex<GparityFlavourIndex>(src_rotated,0,fcol); //only need upper flavor component as X-conjugate
      mixedPrecInvertGen(prop_fcol[fcol], midprop_fcol[fcol], src_fcol, xconj_action_d, xconj_action_f, tol, inner_tol, do_midprop, evals, evecs);
    }    

    get2fXconjMatrix(prop, prop_fcol[0],prop_fcol[1]);
    prop = mulHdagRight(prop);

    if(do_midprop){
      get2fXconjMatrix(midprop, midprop_fcol[0],midprop_fcol[1]);
      midprop = mulHdagRight(midprop);
    }
  }


  //No midprop, with evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  LatticeSCFmatrixD mixedPrecInvertXconj(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, double tol, double inner_tol,
					 std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    LatticeSCFmatrixD tmp(msrc.Grid()), prop(msrc.Grid());
    mixedPrecInvertGenXconj(prop,tmp,msrc,action_d,action_f,tol,inner_tol,false,evals,evecs);
    return prop;
  }
  //No midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  LatticeSCFmatrixD mixedPrecInvertXconj(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, double tol, double inner_tol){
    return mixedPrecInvertXconj(msrc,action_d,action_f,tol,inner_tol,(std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }


  //With midprop and evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       double tol, double inner_tol,
				       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    mixedPrecInvertGenXconj(prop,midprop,msrc,action_d,action_f,tol,inner_tol,true,evals,evecs);
  }    
  //With midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       double tol, double inner_tol){
    mixedPrecInvertWithMidPropXconj(prop, midprop, msrc, action_d, action_f, tol, inner_tol, (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }





};
