#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //Chiral condensate for G-parity: 1/2 ( psibar psi ) = 1/2 ( dbar d + ubar u )
  // 1/2 < psibar psi > = -1/2 tr G  
  // Compute on every site and volume, spin and color avg:
  // -1/(24 V) \sum_x \sum_i G_ii(x,x)
  // Obtain G(x,x) from stochastic random volume source:  
  //     < eta_i(x) eta*_j(y) >_hit ~ delta_ij delta_xy
  // chi_i(x) = \sum_z G_ik(x,z) eta_k(z)
  //     < chi_i(x) eta*_j(y) >_hit ~      < \sum_z G_ik(x,z) eta_k(z)eta*_j(y) >_hit ~  G_ij(x,y)
  // Thus
  // \sum_x \sum_i < chi_i(x) eta*_i(x) >_hit ~ \sum_x \sum_i G_ii(x,x)
  //   = innerProduct( chi*, eta ) *          but we only care about the real part
  RealD chiralCondensate(const FermionFieldD &sol_rnd_vol, const FermionFieldD &src_rnd_vol){
    GridBase* UGrid = sol_rnd_vol.Grid();
    RealD vol = 1.;
    for(int i=0;i<4;i++) vol *= UGrid->GlobalDimensions()[i]; //source is 4d!

    RealD v = innerProduct(conjugate(sol_rnd_vol), src_rnd_vol).real();
    v = v * RealD(1./24/vol); // note I removed the - sign to match CPS conventions, but did not include the F_CLASS_DWF 5-M5 normalization
    return v; 
  }
};
