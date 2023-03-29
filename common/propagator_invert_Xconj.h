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

    T out = mulPplusLeft(T(h1pi*in + h1mi*(sigma3*in))) + mulPminusLeft(T(h1mi*(sigma1*in) - h1mi*(sigma2*in)));
    resetDeviceStackMemory();
    return out;
  }
  template<typename T>
  T mulHdagLeft(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    T out = mulPplusLeft(T(h1mi*in + h1pi*(sigma3*in))) + mulPminusLeft(T(h1pi*(sigma1*in) - h1pi*(sigma2*in)));
    resetDeviceStackMemory();
    return out;
  }
  template<typename T>
  T mulHRight(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    T out = mulPplusRight(T(h1pi*in + h1mi*in*sigma3)) + mulPminusRight(T(h1mi*in*sigma1 - h1mi*in*sigma2));
    resetDeviceStackMemory();
    return out;
  }
  template<typename T>
  inline T mulHdagRight(const T &in){
    static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
    static GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);
    static ComplexD h1pi(0.5,0.5);
    static ComplexD h1mi(0.5,-0.5);

    T out = mulPplusRight(T(h1mi*in + h1pi*in*sigma3)) + mulPminusRight(T(h1pi*in*sigma1 - h1pi*in*sigma2));
    resetDeviceStackMemory();
    return out;
  }




  //Right multiply by U
  template<typename T>
    T mulURight(const T &in){
    static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);

    T out = 0.5*(in + in*Xmatrix()) + 0.5*(in*sigma3 - (in*Xmatrix())*sigma3);
    resetDeviceStackMemory();
    return out;
  }
  
  //Invert using the X-conjugate action
  //NOTE: the source *must* have flavor structure  
  //|  A       B   |
  //| XB*X   -XA*X |
  //which includes
  // | phi   0   |
  // |  0   phi* |  
  //where phi is a spin-color matrix that commutes with X=C g5  
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertGenXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
			       const LatticeSCFmatrixD &msrc, FermionActionD &xconj_action_d, FermionActionF &xconj_action_f, 
			       const MixedCGargs &args,
			       bool do_midprop,
			       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    LatticeSCFmatrixD src_rotated = mulHRight(msrc); //the columns of this are X-conjugate
    LatticePropagatorD src_fcol(msrc.Grid()), tmp(msrc.Grid());
    std::vector<LatticePropagatorD> prop_fcol(2, msrc.Grid()), midprop_fcol(2, msrc.Grid());

    for(int fcol=0;fcol<Ngp;fcol++){
      src_fcol = PeekIndex<GparityFlavourIndex>(src_rotated,0,fcol); //only need upper flavor component as X-conjugate      
      {
	//Check the column is actually X-conjugate!
	tmp = PeekIndex<GparityFlavourIndex>(src_rotated,1,fcol) + Xmatrix()*conjugate(src_fcol);
	assert(norm2(tmp) < 1e-12);
      }
      mixedPrecInvertGen(prop_fcol[fcol], midprop_fcol[fcol], src_fcol, xconj_action_d, xconj_action_f, args, do_midprop, evals, evecs);
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
  LatticeSCFmatrixD mixedPrecInvertXconj(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, const MixedCGargs &args,
					 std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    LatticeSCFmatrixD tmp(msrc.Grid()), prop(msrc.Grid());
    mixedPrecInvertGenXconj(prop,tmp,msrc,action_d,action_f,args,false,evals,evecs);
    return prop;
  }
  //No midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  LatticeSCFmatrixD mixedPrecInvertXconj(const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, const MixedCGargs &args){
    return mixedPrecInvertXconj(msrc,action_d,action_f,args,(std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }


  //With midprop and evecs
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void mixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       const MixedCGargs &args,
				       std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    mixedPrecInvertGenXconj(prop,midprop,msrc,action_d,action_f,args,true,evals,evecs);
  }    
  //With midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void mixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &prop, LatticeSCFmatrixD &midprop, 
				       const LatticeSCFmatrixD &msrc, FermionActionD &action_d, FermionActionF &action_f, 
				       const MixedCGargs &args){
    mixedPrecInvertWithMidPropXconj(prop, midprop, msrc, action_d, action_f, args, (std::vector<Real> const*)nullptr, (std::vector<FermionField1fD> const *)nullptr);
  }


  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertGenXconj(std::vector<LatticeSCFmatrixD> &prop, std::vector<LatticeSCFmatrixD> &midprop, 
					const std::vector<LatticeSCFmatrixD> &msrc, 
					FermionActionD &xconj_action_d, 
					FermionActionD &xconj_subgrid_action_d, FermionActionF &xconj_subgrid_action_f, 
					const MixedCGargs &args,
					bool do_midprop,
					std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    int nsrc = msrc.size();
    GridBase* grid = msrc[0].Grid();
    LatticePropagatorD tmp(grid);
    LatticeSCFmatrixD src_rotated(grid);
    
    int ninv = Ngp * nsrc; //mapping   f + Ngp * s
    std::vector<LatticePropagatorD> src_fcol(ninv,grid);

    for(int s=0;s<nsrc;s++){
      src_rotated = mulHRight(msrc[s]); //the columns of this are X-conjugate

      for(int fcol=0;fcol<Ngp;fcol++){
	int idx = fcol + Ngp * s;
	src_fcol[idx] = PeekIndex<GparityFlavourIndex>(src_rotated,0,fcol); //only need upper flavor component as X-conjugate      
	tmp = PeekIndex<GparityFlavourIndex>(src_rotated,1,fcol) + Xmatrix()*conjugate(src_fcol[idx]);       //Check the columns are actually X-conjugate!
	assert(norm2(tmp) < 1e-12);
      }
    }      

    std::vector<LatticePropagatorD> prop_fcol(ninv, grid), midprop_fcol(do_midprop ? ninv : 0, grid);    
    splitGridMixedPrecInvertGen(prop_fcol, midprop_fcol, src_fcol, xconj_action_d, xconj_subgrid_action_d, xconj_subgrid_action_f, args, do_midprop, evals, evecs);

    for(int s=0;s<nsrc;s++){
      get2fXconjMatrix(prop[s], prop_fcol[0+Ngp*s],prop_fcol[1+Ngp*s]);
      prop[s] = mulHdagRight(prop[s]);

      if(do_midprop){
	get2fXconjMatrix(midprop[s], midprop_fcol[0+Ngp*s],midprop_fcol[1+Ngp*s]);
	midprop[s] = mulHdagRight(midprop[s]);
      }
    }
  }


  //Multi-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertXconj(std::vector<LatticeSCFmatrixD> &msol,
				     const std::vector<LatticeSCFmatrixD> &msrc,
				     FermionActionD &action_d,
				     FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				     const MixedCGargs &args,
				     std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<LatticeSCFmatrixD> tmp;
    splitGridMixedPrecInvertGenXconj(msol, tmp, msrc, action_d, subgrid_action_d, subgrid_action_f, args, false, evals, evecs);
  }
  //Multi-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertXconj(std::vector<LatticeSCFmatrixD> &msol,
				     const std::vector<LatticeSCFmatrixD> &msrc,
				     FermionActionD &action_d,
				     FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				     const MixedCGargs &args){
    splitGridMixedPrecInvertXconj(msol, msrc, action_d, subgrid_action_d, subgrid_action_f, args, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //Single-src, no midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertXconj(LatticeSCFmatrixD &msol,
				     const LatticeSCFmatrixD &msrc,
				     FermionActionD &action_d,
				     FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				     const MixedCGargs &args,
				     std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    
    std::vector<LatticeSCFmatrixD> tmp_msrc(1, msrc);
    std::vector<LatticeSCFmatrixD> tmp_msol;
    std::vector<LatticeSCFmatrixD> tmp_msol_mid;
    splitGridMixedPrecInvertGenXconj(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, args, false, evals, evecs);
    msol = tmp_msol[0];
  }
  //Single-src, no midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertXconj(LatticeSCFmatrixD &msol,
				     const LatticeSCFmatrixD &msrc,
				     FermionActionD &action_d,
				     FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
				     const MixedCGargs &args){
    splitGridMixedPrecInvertXconj(msol,msrc,action_d,subgrid_action_d,subgrid_action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }
				


  //Multi-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidPropXconj(std::vector<LatticeSCFmatrixD> &msol, std::vector<LatticeSCFmatrixD> &msol_mid,
						const std::vector<LatticeSCFmatrixD> &msrc,
						FermionActionD &action_d,
						FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
						const MixedCGargs &args,
						std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    splitGridMixedPrecInvertGenXconj(msol, msol_mid, msrc, action_d, subgrid_action_d, subgrid_action_f, args, true, evals, evecs);
  }
  //Multi-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidPropXconj(std::vector<LatticeSCFmatrixD> &msol, std::vector<LatticeSCFmatrixD> &msol_mid,
						const std::vector<LatticeSCFmatrixD> &msrc,
						FermionActionD &action_d,
						FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
						const MixedCGargs &args){
    splitGridMixedPrecInvertWithMidPropXconj(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,args,(std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


  //Single-src with midprop
  template<typename FermionActionD, typename FermionActionF, typename EvecFieldType>
  void splitGridMixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						const LatticeSCFmatrixD &msrc,
						FermionActionD &action_d,
						FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
						const MixedCGargs &args,
						std::vector<Real> const* evals, std::vector<EvecFieldType> const * evecs){
    std::vector<LatticeSCFmatrixD> tmp_msrc(1, msrc);
    std::vector<LatticeSCFmatrixD> tmp_msol;
    std::vector<LatticeSCFmatrixD> tmp_msol_mid;        
    splitGridMixedPrecInvertGenXconj(tmp_msol, tmp_msol_mid, tmp_msrc, action_d, subgrid_action_d, subgrid_action_f, args, true, evals, evecs);
    msol = tmp_msol[0];
    msol_mid = tmp_msol_mid[0];
  }
  //Single-src with midprop, no evecs
  template<typename FermionActionD, typename FermionActionF>
  void splitGridMixedPrecInvertWithMidPropXconj(LatticeSCFmatrixD &msol, LatticeSCFmatrixD &msol_mid,
						const LatticeSCFmatrixD &msrc,
						FermionActionD &action_d,
						FermionActionD &subgrid_action_d, FermionActionF &subgrid_action_f,
						const MixedCGargs &args){
    splitGridMixedPrecInvertWithMidPropXconj(msol,msol_mid,msrc,action_d,subgrid_action_d,subgrid_action_f,args, (std::vector<Real> const*)nullptr, (std::vector<typename FermionActionD::FermionField> const *)nullptr);
  }


};
