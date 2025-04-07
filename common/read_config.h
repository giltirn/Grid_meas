#pragma once

#include<Grid/Grid.h>

namespace GridMeas{
  using namespace Grid;

  template<typename GimplD>
  void readConfiguration(LatticeGaugeFieldD &U,
			 GridSerialRNG &sRNG,
			 GridParallelRNG &pRNG,
			 int traj, 
			 const std::string &cfg_stub,
			 const std::string &rng_stub){

    std::string dummy;
    CheckpointerParameters p(cfg_stub, dummy, rng_stub);
    NerscHmcCheckpointer<GimplD> cp(p);

    cp.CheckpointRestore(traj, U, sRNG, pRNG);
  }

  //For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
  template<typename GimplD>
  void readCPSconfiguration(LatticeGaugeFieldD &U,
			    GridSerialRNG &sRNG,
			    GridParallelRNG &pRNG,
			    int traj, 
			    const std::string &cfg_stub){

    NerscIO::exitOnReadPlaquetteMismatch() = false;

    CheckpointerParameters p(cfg_stub, "pooh");
    NerscHmcCheckpointer<GimplD> cp(p);
    typedef GaugeStatistics<GimplD> GaugeStats;
    
    std::string config, dummy;
    cp.build_filenames(traj, p, config, dummy, dummy);
       
    cp.check_filename(config);
  
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, config);

    NerscIO::exitOnReadPlaquetteMismatch() = true;

    std::vector<int> seeds4({traj, traj+1, traj+2, traj+3});
    pRNG.SeedFixedIntegers(seeds4);
    sRNG.SeedFixedIntegers(seeds4); 
  }

  template<typename GimplD>
  void readConfiguration(LatticeGaugeFieldD &U,
			 const std::string &filename){
    typedef GaugeStatistics<GimplD> GaugeStats;    
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, filename);
  }

  //For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
  template<typename GimplD>
  void readCPSconfiguration(LatticeGaugeFieldD &U,
			    const std::string &filename){

    NerscIO::exitOnReadPlaquetteMismatch() = false;
    typedef GaugeStatistics<GimplD> GaugeStats;
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, filename);
    NerscIO::exitOnReadPlaquetteMismatch() = true;
  }


};
