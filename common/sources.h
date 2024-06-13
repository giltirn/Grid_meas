#pragma once

#include "defines.h"
#include "utils.h"
#include "momentum.h"

namespace GridMeas{
  using namespace Grid;

  template<typename Action, typename std::enable_if<std::is_same<typename GridTypeMapper<typename Action::SitePropagator>::scalar_objectD , SCFmatrixD>::value, int>::type = 0>
  SCFmatrixD unitSiteSrc(){
    SCFmatrixD out;
    out = Zero();
    for(int f=0;f<Ngp;f++)
      for(int s=0;s<Ns;s++)  
	for(int c=0;c<Nc;c++)  
	  out(f,f)(s,s)(c,c) = 1.;
    return out;
  }
  template<typename Action, typename std::enable_if<std::is_same<typename GridTypeMapper<typename Action::SitePropagator>::scalar_objectD, SpinColourMatrixD>::value, int>::type = 0>
  SpinColourMatrixD unitSiteSrc(){
    SpinColourMatrixD out;
    out = Zero();
    for(int s=0;s<Ns;s++)  
      for(int c=0;c<Nc;c++)  
	out()(s,s)(c,c) = 1.;
    return out;
  }


  template<typename ActionD>
  typename ActionD::PropagatorField pointSource(const Coordinate &coord, GridBase* UGrid){
    auto one = unitSiteSrc<ActionD>();
  
    typename ActionD::PropagatorField src(UGrid);
    src = Zero();
    pokeSite(one, src, coord);
    return src;
  }
  template<typename ActionD>
  typename ActionD::PropagatorField Z2wallSource(const int t, GridParallelRNG &RNG, GridBase* UGrid){
    auto one = unitSiteSrc<ActionD>();
    decltype(one) mone = -one;

    typename ActionD::PropagatorField zone(UGrid), zmone(UGrid), zero(UGrid);
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
  
    typename ActionD::PropagatorField src = where(rfield < RealD(0.5), zmone, zone);  
    return src;
  }
  template<typename ActionD>
  typename ActionD::PropagatorField wallSource(const int t, GridBase* UGrid){
    auto one = unitSiteSrc<ActionD>();

    typename ActionD::PropagatorField zone(UGrid), zero(UGrid);
    zone = one;
    zero = Zero();
  
    //Constrain to wall
    LatticeInteger coord_t(UGrid);
    LatticeCoordinate(coord_t, 3);
    zone = where(coord_t == Integer(t), zone, zero);

    return zone;
  }

  //\eta(x,tau) = \delta_{tau,t} exp(-i \vec p \cdot \vec x)
  template<typename ActionD>
  typename ActionD::PropagatorField momentumWallSource(const std::vector<double> &p, const int t, GridBase* UGrid){    
    LatticeComplexD phase_field = phaseField(p, UGrid); //exp(-i \vec p \cdot \vec x)
    typename ActionD::PropagatorField src = wallSource<ActionD>(t,UGrid);
    src = phase_field * src;
    return src;
  }

  //\eta(x,tau) = \delta_{tau,t} \prod_{i=0}^3 cos(p_i x_i)
  template<typename ActionD>
  typename ActionD::PropagatorField cosineWallSource(const std::vector<double> &p, const int t, GridBase* UGrid){    
    LatticeComplexD phase_field = cosineMomentumField(p, UGrid);
    typename ActionD::PropagatorField src = wallSource<ActionD>(t,UGrid);
    src = phase_field * src;
    return src;
  }

  //G-parity only
  FermionFieldD randomGaussianVolumeSource(GridParallelRNG &rng, GridBase* UGrid){
    FermionFieldD out(UGrid);
    //Default sigma^2 = 1     CPS uses sigma^2=1/2 for random source (cf alg_pbp.C)
    //The reason for choosing sigma^2=1/2 is because we want 1=<x* x> = <x_r^2> + <x_i^2>     and for mu=0  <x_r^2> = <x_i>^2 = sigma^2 

    //Properties of gaussian distribution   if   X ~ N(0,1)   then   aX  ~ N(0, |a|)
    gaussian(rng, out);
    out = out * RealD(sqrt(0.5)); //re,im
    return out;
  }

  //\eta(x,tau) = \delta_{tau,t} [ P_p e^{-i\vec p\cdot\vec x} + (P_p)^*  e^{i\vec p\cdot\vec x} ]
  //where P_p is the projector suitable for a \bar\psi field with momentum p
  //G-parity only
  LatticeSCFmatrixD gparityCosineWallSource(const std::vector<int> &p, const int t, GridBase* UGrid){
    GparityFlavour P_p(getProjector(p, true));
    
    LatticeSCFmatrixD mom_wall = momentumWallSource<CayleyFermion5D<GparityWilsonImplD> >(getPhysicalMomentum(p),t,UGrid);
    LatticeSCFmatrixD out = mom_wall * P_p;
    out = out + conjugate(out);
    return out;
  }
  template<typename ActionD>
  typename ActionD::PropagatorField volumeSource(GridBase* UGrid){
    auto one = unitSiteSrc<ActionD>();

    typename ActionD::PropagatorField zone(UGrid);
    zone = one;
    return zone;
  }
  //\eta(x) = exp(-i p_\mu x_mu) {  I_{24x24} or I_{12x12}  }
  template<typename ActionD>
  typename ActionD::PropagatorField fourMomentumVolumeSource(const std::vector<double> &p, GridBase* UGrid){    
    assert(p.size() == 4);
    LatticeComplexD phase_field = phaseFieldFour(p, UGrid);
    typename ActionD::PropagatorField src = volumeSource<ActionD>(UGrid);
    src = phase_field * src;
    return src;
  }

};
