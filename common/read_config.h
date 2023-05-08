#pragma once

#include<Grid/Grid.h>

namespace GridMeas{
  using namespace Grid;

  
  void readConfiguration(LatticeGaugeFieldD &U,
			 GridSerialRNG &sRNG,
			 GridParallelRNG &pRNG,
			 int traj, 
			 const std::string &cfg_stub,
			 const std::string &rng_stub){

    CheckpointerParameters p(cfg_stub, rng_stub);
    NerscHmcCheckpointer<ConjugateGimplD> cp(p);

    cp.CheckpointRestore(traj, U, sRNG, pRNG);
  }

  //For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
  void readCPSconfiguration(LatticeGaugeFieldD &U,
			    GridSerialRNG &sRNG,
			    GridParallelRNG &pRNG,
			    int traj, 
			    const std::string &cfg_stub){

    NerscIO::exitOnReadPlaquetteMismatch() = false;

    CheckpointerParameters p(cfg_stub, "pooh");
    NerscHmcCheckpointer<ConjugateGimplD> cp(p);
    typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
    
    std::string config, dummy;
    cp.build_filenames(traj, p, config, dummy);  
    cp.check_filename(config);
  
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, config);

    NerscIO::exitOnReadPlaquetteMismatch() = true;

    std::vector<int> seeds4({traj, traj+1, traj+2, traj+3});
    pRNG.SeedFixedIntegers(seeds4);
    sRNG.SeedFixedIntegers(seeds4); 
  }


  void readConfiguration(LatticeGaugeFieldD &U,
			 const std::string &filename){
    typedef GaugeStatistics<ConjugateGimplD> GaugeStats;    
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, filename);
  }

  //For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
  void readCPSconfiguration(LatticeGaugeFieldD &U,
			    const std::string &filename){

    NerscIO::exitOnReadPlaquetteMismatch() = false;
    typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, filename);
    NerscIO::exitOnReadPlaquetteMismatch() = true;
  }


};
