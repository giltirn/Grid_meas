#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;
  
  //Compute t^2 E(t) for 0<t<=Nstep*epsilon in increments of epsilon
  //Return vector of  [ t, t^2 E(t) ]
  //If V != nullptr, the smeared field is placed in V
  std::vector<std::pair<RealD,RealD> > WilsonFlowEnergyDensity(const int Nstep, const double epsilon, const LatticeGaugeFieldD &U, LatticeGaugeFieldD *V = nullptr){
    WilsonFlow<ConjugateGimplD> wflow(Nstep, epsilon);
    std::vector<RealD> vals = V != nullptr ? wflow.flowMeasureEnergyDensityCloverleaf(*V,U) : wflow.flowMeasureEnergyDensityCloverleaf(U);
    assert(vals.size() == Nstep );
    std::vector<std::pair<RealD, RealD> > out(Nstep);
    for(int i=0;i<Nstep;i++){
      out[i] = { epsilon*(i+1), vals[i] };
    }
    return out;
  }
};
