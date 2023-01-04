#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  RealD avgOrientedPlaquette(const int mu, const int nu, const LatticeGaugeFieldD &Umu){
    std::vector<LatticeColourMatrixD> U(Nd, Umu.Grid());
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }
    double vol = Umu.Grid()->gSites();

    LatticeComplexD Plaq(Umu.Grid());
    WilsonLoops<ConjugateGimplD>::traceDirPlaquette(Plaq, U, mu, nu);
    auto Tp = sum(Plaq);
    auto p = TensorRemove(Tp);
    RealD pr = p.real();
    return pr / vol / Nc;
  }

};


