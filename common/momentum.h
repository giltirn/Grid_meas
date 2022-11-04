#pragma once

#include<Grid/Grid.h>

namespace GridMeas{
  using namespace Grid;

  //f(\vec x, t) = exp(-i \vec p \cdot \vec x)
  LatticeComplexD phaseField(const std::vector<double> &p, GridBase* Grid){
    LatticeComplexD out(Grid);
    LatticeComplexD coor(Grid);
    out=Zero();
    for(int mu=0;mu<3;mu++){
      LatticeCoordinate(coor,mu);
      out = out + p[mu] * coor;
    }
    ComplexD ci(0.0,-1.0);
    out = exp(out*ci);
    return out;
  }

  //Phase field for a cosine source, cf https://rbc.phys.columbia.edu/rbc_ukqcd/phd_thesis/Lightman%20Matthew/thesis_v2_0_distribute.pdf  eq 6.18
  //f(\vec x, t) = \prod_{mu=0}^{2} cos(p_\mu x_\mu)
  LatticeComplexD cosineMomentumField(const std::vector<double> &p, GridBase* Grid){
    LatticeComplexD out(Grid);
    LatticeComplexD coor(Grid);

    LatticeCoordinate(coor,0);
    out=cos(p[0]*coor);
    for(int mu=1;mu<3;mu++){
      LatticeCoordinate(coor,mu);
      out = out * cos(p[mu] * coor);
    }
    return out;
  }

  //exp(-i \vec p \cdot x0)
  ComplexD phase(const std::vector<double> &p, const Coordinate &x0){
    ComplexD out(0);
    for(int mu=0;mu<3;mu++){
      out = out + p[mu] * x0[mu];
    }
    ComplexD ci(0.0,-1.0);
    out = exp(out*ci);
    return out;
  }

  //Get the physical momentum in lattice units. Input units are :  2pi/L  (non-Gparity)   pi/2L (G-parity)
  std::vector<double> getPhysicalMomentum(const std::vector<int> &p){
    assert(p.size() == 3);
    std::vector<int> gpdirs = ConjugateGimplD::getDirections();
    Coordinate latt = GridDefaultLatt();
    std::vector<double> pphys(3,0);
    for(int i=0;i<3;i++){
      double punit = gpdirs[i] ? M_PI/2./latt[i] : 2*M_PI/latt[i];
      pphys[i] = p[i]*punit;
    }
    return pphys;
  }


  //Get the projector 1/2(1\pm\sigma_2) appropriate for a given quark momentum
  //p should be in units of pi/2L
  GparityFlavour::Algebra getProjector(const std::vector<int> &p){
    //psi_+  :   p = (..., -7, -3, 1, 5, 9, ...) pi/2L
    //psi_-  :   p = (..., -9, -5, -1, 3, 7, ...) pi/2L

    //are either described as (1+4n) or -(1+4n) for integer n
    std::vector<int> gpdirs = ConjugateGimplD::getDirections();
    int pm = -99;
    for(int mu=0;mu<3; mu++){
      if(gpdirs[mu]){
	assert( abs(p[mu]) % 2 == 1 ); //odd
	int pm_mu;
	//try (1+4n)
	int n = (p[mu] - 1) / 4;
	if(4*n + 1 == p[mu]){
	  pm_mu = 1;
	
	}else{
	  //try -(1+4n)
	  n = (-p[mu] - 1) / 4;
	  assert( -(4*n + 1) == p[mu] );
	  pm_mu = -1;
	}
    
	if(pm == -99) pm = pm_mu;
	else assert( pm == pm_mu ); //must be consistent
      }
    }
    assert(pm != -99);
  
    return pm == 1 ? GparityFlavour::Algebra::ProjPlus : GparityFlavour::Algebra::ProjMinus;
  }


  std::string momstr(const std::vector<int> &p){
    std::stringstream ss;
    for(int i=0;i<3;i++){
      if(p[i] < 0) ss << "_" << (-p[i]);
      else ss << p[i];
    }
    return ss.str();
  }


};
