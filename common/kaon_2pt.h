#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //https://arxiv.org/pdf/1908.08640.pdf eq 179
  //momentum wall sources

  //C(t) =   1/2 \sum_{\vec x, \vec y_1, \vec y_2} e^{i( \vec p\cdot \vec x - \vec p_1 \cdot \vec y_1 - \vec p_2 \cdot \vec y_2 ) } tr [ \gamma^5 G_l(\vec x,t; \vec y_1,t_0)\eta(\vec y_1,t_0) \gamma^5 \Xi  \eta(\vec y_2,t_0)G_h(\vec y_2,t_0; \vec x, t) ]
  //\vec p_1 + \vec p_2 = \vec p = 0

  //\Xi is a projector 1/2(1+-sigma_2)

  //Use \vec p_2 = -\vec p_1

  //C(t) =   1/2 \sum_{\vec x} tr [ \gamma^5 R_l(\vec x,t)\eta(\vec y_1,t_0) \gamma^5 \Xi  S_h(\vec x, t) ]

  //R_l(\vec x, t; \vec p_1) = \sum_{\vec y_1} e^{-i \vec p_1 \cdot \vec y_1} G_l(\vec x, t; \vec y_1, t_0)\eta(\vec y_1, t_0)
  //S_h(\vec x, t; \vec p_1) = \sum_{\vec y_2} e^{+i \vec p_1 \cdot \vec y_2} \eta(\vec y_2, t_0) G_h(\vec y_2, t_0 ; \vec x, t )

  //Iff \eta(\vec y_2, t_0) has diagonal spin, color, flavor structure and *is real*
  //S_h(\vec x, t; \vec p_1) = \gamma^5 [ \sum_{\vec y_2} e^{-i \vec p_1 \cdot \vec y_2}   G_h(\vec x, t; \vec y_2, t_0 ) \eta(\vec y_2, t_0) ]^\dagger \gamma^5
  //                         = \gamma^5 [ R_h(\vec x, t; \vec p_1) ]^\dagger \gamma^5
  
  //Where R_h(\vec x, t; \vec p_1) = \sum_{\vec y_2} e^{-i \vec p_1 \cdot \vec y_2} G_h(\vec x, t; \vec y_2, t_0)\eta(\vec y_2, t_0)

  //C(t) =   1/2 \sum_{\vec x} tr [ R_l(\vec x,t; \vec p_1) \Xi  R_h(\vec x,t; \vec p_1)^\dagger ]

  //\Xi should be the projector for momentum \vec p_2 = -\vec p_1

  std::vector<RealD> momWallSourceKaonCorrelator(const std::vector<int> &p1, int t0, const LatticeSCFmatrixD& R_l_p1, const LatticeSCFmatrixD& R_h_p1){
    std::cout << GridLogMessage << "Starting momentum wall source kaon correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R_l_p1.Grid();
  
    std::vector<int> p2 = {-p1[0], -p1[1], -p1[2]};
    GparityFlavour Xi(getProjector(p2));
    std::cout << GridLogMessage << "Projector for p2=" << p2 << " : " << getProjector(p2) << std::endl;

    LatticeSCFmatrixD sqb = R_l_p1 * Xi * adj(R_h_p1);
  
    std::vector<SCFmatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    for(int t=0;t<Lt;t++)
      std::cout << GridLogMessage << "C(t)=" << 0.5*trace(Ctm[t]) << std::endl;
    
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = 0.5 * real( trace(Ctm[t]) );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }

};
