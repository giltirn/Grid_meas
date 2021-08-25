#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;
  
  inline void LandauGaugeFix(LatticeGaugeFieldD &U, double alpha=0.1, double stop=1e-10, bool fourier_accelerate=true){
    std::cout << GridLogMessage << "Gauge fixing to Landau gauge" << std::endl;
    FourierAcceleratedGaugeFixer<ConjugateGimplD>::SteepestDescentGaugeFix(U,alpha,100000,stop, stop,fourier_accelerate);
  }
  inline void CoulombGaugeFix(LatticeGaugeFieldD &U, double alpha=0.1, double stop=1e-10, bool fourier_accelerate=true){
    std::cout << GridLogMessage << "Gauge fixing to Coulomb gauge" << std::endl;
    FourierAcceleratedGaugeFixer<ConjugateGimplD>::SteepestDescentGaugeFix(U,alpha,100000,stop, stop,fourier_accelerate,Nd-1);
  }

  //Write the gauge fixed configuration to a file  filename_stub + ".traj"
  inline void writeGaugeFixedConfiguration(const LatticeGaugeFieldD &U, const std::string &filename_stub, int traj){
    std::string file = filename_stub + "." + std::to_string(traj);
    std::cout << GridLogMessage << "Writing gauge fixed configuration to " << file << std::endl;
    typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
    int precision32 = 1;
    int tworow = 0;
    NerscIO::writeConfiguration<GaugeStats>(const_cast<LatticeGaugeFieldD &>(U), file, tworow, precision32);
  }

  //Write the gauge fixed configuration from a file  filename_stub + ".traj"
  inline void readGaugeFixedConfiguration(LatticeGaugeFieldD &U, const std::string &filename_stub, int traj){
    std::string file = filename_stub + "." + std::to_string(traj);
    std::cout << GridLogMessage << "Reading gauge fixed configuration from " << file << std::endl;
    typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
    FieldMetaData header;
    NerscIO::readConfiguration<GaugeStats>(U, header, file);
  }
 

};
