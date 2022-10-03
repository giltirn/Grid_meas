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
  //       topq: [ t, [Q(tslice=0), Q(tslice=1)... Q(tslice=Lt-1)] ]
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

  struct WilsonFlowIO{
    bool do_energy_density_clover;
    std::vector<std::pair<RealD,RealD> > energy_density_clover; //[ t, t^2 E(t) ]

    bool do_energy_density_plaq;
    std::vector<std::pair<RealD,RealD> > energy_density_plaq; //[ t, t^2 E(t) ]
    
    bool do_timeslice_topq;
    int timeslice_topq_meas_freq;
    std::vector<std::pair<RealD,std::vector<RealD> > > timeslice_topq; //[ t, [Q(tslice=0), Q(tslice=1)... Q(tlice=Lt-1)] ]

    bool do_timeslice_plaq;
    int timeslice_plaq_meas_freq;
    std::vector<std::pair<RealD,std::vector<RealD> > > timeslice_plaq; //[ t, [P(tslice=0), P(tslice=1)... P(tslice=Lt-1)] ]

    WilsonFlowIO(){
      do_energy_density_clover = true;
      do_energy_density_plaq = false;
      do_timeslice_topq = false;
      do_timeslice_plaq = false;
      timeslice_topq_meas_freq=1;
      timeslice_plaq_meas_freq=1;
    }
    //Clear data outputs
    void clear(){
      energy_density_clover.clear();
      energy_density_plaq.clear();
      timeslice_topq.clear();
      timeslice_plaq.clear();
    }

  };
  
  //Wilson flow smearing with configurable outputs
  void WilsonFlowMeasGeneralOpt(WilsonFlowIO &meas,
				const int Nstep, const double epsilon, bool do_adaptive, const double adaptive_maxTau, const double adaptive_tol,
				const LatticeGaugeFieldD &U, LatticeGaugeFieldD *V = nullptr){
    WilsonFlowBase<ConjugateGimplD>* wflow = nullptr;
    if(do_adaptive) wflow = new WilsonFlowAdaptive<ConjugateGimplD>(epsilon, adaptive_maxTau, adaptive_tol);
    else wflow = new WilsonFlow<ConjugateGimplD>(epsilon, Nstep);

    wflow->resetActions();    
    meas.clear();

    if(meas.do_energy_density_clover)
      wflow->addMeasurement(1, [&wflow,&meas](int step, RealD t, const LatticeGaugeField &U){ 
	  std::cout << GridLogMessage << "[WilsonFlow] Computing Cloverleaf energy density for step " << step << std::endl;
	  meas.energy_density_clover.push_back( {t, wflow->energyDensityCloverleaf(t,U)} );
	});

    if(meas.do_energy_density_plaq)
      wflow->addMeasurement(1, [&wflow,&meas](int step, RealD t, const LatticeGaugeField &U){ 
	  std::cout << GridLogMessage << "[WilsonFlow] Computing Plaquette energy density for step " << step << std::endl;
	  meas.energy_density_plaq.push_back( {t, wflow->energyDensityPlaquette(t,U)} );
	});
    
    if(meas.do_timeslice_topq)
      wflow->addMeasurement(meas.timeslice_topq_meas_freq, [&meas](int step, RealD t, const LatticeGaugeField &U){ 
	  std::cout << GridLogMessage << "[WilsonFlow] Computing timeslice topologial charge for step " << step << std::endl;
	  meas.timeslice_topq.push_back( {t, WilsonLoops<ConjugateGimplD>::TimesliceTopologicalCharge5Li(U) } );
	});      

    if(meas.do_timeslice_plaq)
      wflow->addMeasurement(meas.timeslice_plaq_meas_freq, [&meas](int step, RealD t, const LatticeGaugeField &U){ 
	  std::cout << GridLogMessage << "[WilsonFlow] Computing timeslice plaquette charge for step " << step << std::endl;
	  meas.timeslice_plaq.push_back( {t, WilsonLoops<ConjugateGimplD>::timesliceAvgSpatialPlaquette(U) } );
	});      

    LatticeGaugeFieldD Vtmp(U.Grid());
    wflow->smear(Vtmp, U);
    delete wflow;

    if(V) *V = std::move(Vtmp);
  }
  inline void WilsonFlowMeasGeneral(WilsonFlowIO &meas,
			     const int Nstep, const double epsilon,
			     const LatticeGaugeFieldD &U, LatticeGaugeFieldD *V = nullptr){
    WilsonFlowMeasGeneralOpt(meas,Nstep,epsilon,false,0,0,U,V);
  }
  inline void WilsonFlowAdaptiveMeasGeneral(WilsonFlowIO &meas,
					    const double epsilon, const double maxTau, const double tol,
					    const LatticeGaugeFieldD &U, LatticeGaugeFieldD *V = nullptr){
    WilsonFlowMeasGeneralOpt(meas,0,epsilon,true,maxTau,tol,U,V);
  }



};
