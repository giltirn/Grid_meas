#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //Generate the contributions to the 5Li topological charge from Wilson loops of the following sizes
  //Use coefficients from hep-lat/9701012
  //1x1 : c1=(19.-55.*c5)/9.
  //2x2 : c2=(1-64.*c5)/9.
  //1x2 : c3=(-64.+640.*c5)/45.
  //1x3 : c4=1./5.-2.*c5
  //3x3 : c5=1./20.
  //Output array contains the loops in the above order
  //V should be a smeared field
  inline std::vector<Real> topologicalCharge5LiContributions(const LatticeGaugeFieldD &V){
    return WilsonLoops<ConjugateGimplD>::TopologicalCharge5LiContributions(V);
  }
  
  inline Real topologicalCharge5Li(const std::vector<Real> &contribs){
    double c5=1./20.;
    double c4=1./5.-2.*c5;
    double c3=(-64.+640.*c5)/45.;
    double c2=(1-64.*c5)/9.;
    double c1=(19.-55.*c5)/9.;

    double Q = c1*contribs[0] + c2*contribs[1] + c3*contribs[2] + c4*contribs[3] + c5*contribs[4];
    return Q;
  }

  inline Real topologicalCharge5Li(const LatticeGaugeFieldD &V){
    return topologicalCharge5Li( topologicalCharge5LiContributions(V) );
  }
};
