#pragma once

#include "defines.h"
#include "utils.h"

namespace GridMeas{
  using namespace Grid;

  SCFmatrixD unitSiteSrc(){
    SCFmatrixD out;
    out = Zero();
    for(int f=0;f<Ngp;f++)
      for(int s=0;s<Ns;s++)  
	for(int c=0;c<Nc;c++)  
	  out(f,f)(s,s)(c,c) = 1.;
    return out;
  }

  LatticeSCFmatrixD pointSource(const Coordinate &coord, GridBase* UGrid){
    SCFmatrixD one = unitSiteSrc();
  
    LatticeSCFmatrixD src(UGrid);
    src = Zero();
    pokeSite(one, src, coord);
    return src;
  }

  LatticeSCFmatrixD Z2wallSource(const int t, GridParallelRNG &RNG, GridBase* UGrid){
    SCFmatrixD one = unitSiteSrc();
    SCFmatrixD mone = -one;

    LatticeSCFmatrixD zone(UGrid), zmone(UGrid), zero(UGrid);
    zone = one;
    zmone = mone;
    zero = Zero();
  
    //Constrain to wall
    LatticeInteger coord_t(UGrid);
    LatticeCoordinate(coord_t, 3);
    zone = where(coord_t == Integer(t), zone, zero);
    zmone = where(coord_t == Integer(t), zmone, zero);

    //Random choice
    LatticeRealD rfield(UGrid);
    random(RNG, rfield); //  0 <=  f(x) < 1
  
    LatticeSCFmatrixD src = where(rfield < RealD(0.5), zmone, zone);  
    return src;
  }

  LatticeSCFmatrixD wallSource(const int t, GridBase* UGrid){
    SCFmatrixD one = unitSiteSrc();

    LatticeSCFmatrixD zone(UGrid), zero(UGrid);
    zone = one;
    zero = Zero();
  
    //Constrain to wall
    LatticeInteger coord_t(UGrid);
    LatticeCoordinate(coord_t, 3);
    zone = where(coord_t == Integer(t), zone, zero);

    return zone;
  }

  FermionFieldD randomGaussianVolumeSource(GridParallelRNG &rng, GridBase* UGrid){
    FermionFieldD out(UGrid);
    //Default sigma^2 = 1     CPS uses sigma^2=1/2 for random source (cf alg_pbp.C)
    //The reason for choosing sigma^2=1/2 is because we want 1=<x* x> = <x_r^2> + <x_i^2>     and for mu=0  <x_r^2> = <x_i>^2 = sigma^2 

    //Properties of gaussian distribution   if   X ~ N(0,1)   then   aX  ~ N(0, |a|)
    gaussian(rng, out);
    out = out * RealD(sqrt(0.5)); //re,im
    return out;
  }




};
