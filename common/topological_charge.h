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

  //outer index is the loop index, inner is t
  inline std::vector<std::vector<Real> > timesliceTopologicalCharge5LiContributions(const LatticeGaugeFieldD &V){
    return WilsonLoops<ConjugateGimplD>::TimesliceTopologicalCharge5LiContributions(V);
  }

  //contributions from timeslice contributions
  inline std::vector<Real> topologicalCharge5LiContributions(const std::vector<std::vector<Real> > &tslice_contribs){
    std::vector<Real> out(tslice_contribs.size(),0);
    for(int i=0;i<tslice_contribs.size();i++)
      for(int t=0;t<tslice_contribs[i].size();t++)
	out[i] += tslice_contribs[i][t];
    return out;
  }
  
  inline std::vector<Real> timesliceTopologicalCharge5Li(const std::vector<std::vector<Real> > &contribs){
    double c5=1./20.;
    double c4=1./5.-2.*c5;
    double c3=(-64.+640.*c5)/45.;
    double c2=(1-64.*c5)/9.;
    double c1=(19.-55.*c5)/9.;

    int Lt = contribs[0].size();
    std::vector<Real> out(Lt,0.);
    for(int t=0;t<Lt;t++)
      out[t] += c1*contribs[0][t] + c2*contribs[1][t] + c3*contribs[2][t] + c4*contribs[3][t] + c5*contribs[4][t];
    return out;
  }

  inline std::vector<Real> timesliceTopologicalCharge5Li(const LatticeGaugeFieldD &V){
    return timesliceTopologicalCharge5Li( timesliceTopologicalCharge5LiContributions(V) );
  }



};
