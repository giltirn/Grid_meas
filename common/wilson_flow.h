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

  //Compute t^2 E(t) for 0<t<=Nstep*epsilon in increments of epsilon
  //and also compute the timeslice topological charge every 'meas_freq' iterations
  //Output: 
  //       energy_density:  [ t, t^2 E(t) ]
  //       topq: [ t, [Q(t=0), Q(t=1)... Q(t=Lt-1)] ]
  //If V != nullptr, the smeared field is placed in V  
  void WilsonFlowEnergyDensityAndTimesliceTopQ(std::vector<std::pair<RealD,RealD> > &energy_density,
					       std::vector<std::pair<RealD,std::vector<RealD> > > &topq,
					       const int Nstep, const double epsilon, const int meas_freq,
					       const LatticeGaugeFieldD &U, LatticeGaugeFieldD *V = nullptr){
    WilsonFlow<ConjugateGimplD> wflow(Nstep, epsilon);
    wflow.resetActions();    
    energy_density.resize(0);
    topq.resize(0);

    wflow.addMeasurement(1, [&wflow,&energy_density](int step, RealD t, const LatticeGaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing Cloverleaf energy density for step " << step << std::endl;
      energy_density.push_back( {t, wflow.energyDensityCloverleaf(t,U)} );
    });      

    wflow.addMeasurement(meas_freq, [&topq](int step, RealD t, const LatticeGaugeField &U){ 
      std::cout << GridLogMessage << "[WilsonFlow] Computing timeslice topologial charge for step " << step << std::endl;
      topq.push_back( {t, WilsonLoops<ConjugateGimplD>::TimesliceTopologicalCharge5Li(U) } );
    });      

    LatticeGaugeFieldD Vtmp(U.Grid());
    wflow.smear(Vtmp, U);
    if(V) *V = std::move(Vtmp);
  }


};
