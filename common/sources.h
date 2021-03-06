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




};
