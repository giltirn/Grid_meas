#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //C(t) =   \sum_{\vec x, \vec y_1, \vec y_2} e^{-i( \vec p\cdot \vec x - \vec p_1 \cdot \vec y_1 - \vec p_2 \cdot \vec y_2 ) } tr [ \gamma^5\sigma_3 G(\vec x, \vec y_1) \gamma^5\sigma_3 \phi(\vec y_1, \vec y_2)  G(\vec y_2, \vec x) ]
  //\vec p_1 + \vec p_2 = \vec p

  //For point sources:
  //C(t) =   \sum_{\vec x} e^{-i\vec p\cdot( \vec x - \vec y_0 ) } tr [ \gamma^5\sigma_3 G(\vec x, \vec y_0) \gamma^5\sigma_3 G(\vec y_0, \vec x) ]
  //     =   \sum_{\vec x} e^{-i\vec p\cdot( \vec x - \vec y_0 ) } tr [ \gamma^5\sigma_3 G(\vec x, \vec y_0) \gamma^5\sigma_3 \gamma^5 G^\dagger(\vec x, \vec y_0) \gamma^5 ]
  //     =   \sum_{\vec x} e^{-i\vec p\cdot \vec x}e^{i\vec p\cdot \vec y_0} tr [ \sigma_3 G(\vec x, \vec y_0) \sigma_3 G^\dagger(\vec x, \vec y_0) ]

  //No need for source flavor projection as not a smeared source

  std::vector<RealD> pointSourcePionCorrelator(const std::vector<int> &p, const Coordinate &y0, const LatticeSCFmatrixD& G){
    std::vector<double> pphys = getPhysicalMomentum(p);
    ComplexD src_phase = conj( phase(pphys,y0 ) ); //exp(i \vec p \cdot y_0 )
    LatticeComplexD snk_phase_field = phaseField(pphys, G.Grid()); //exp(-i \vec p \cdot \vec x)
  
    GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
    Gamma gamma5(Gamma::Algebra::Gamma5);

    LatticeSCFmatrixD sqb =  snk_phase_field *(sigma3* G * sigma3 * adj(G) );
  
    std::vector<SCFmatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace(Ctm[t]) * src_phase );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - y0[3] + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }


  //momentum wall sources
  //https://arxiv.org/pdf/1908.08640.pdf  eq 145 but with pseudoscalar not axial sink
  //C(t) =   1/2 \sum_{\vec x, \vec y_1, \vec y_2} e^{i( \vec p\cdot \vec x - \vec p_1 \cdot \vec y_1 - \vec p_2 \cdot \vec y_2 ) } tr [ \gamma^5\sigma_3 G(\vec x,t; \vec y_1,t_0)\eta(\vec y_1,t_0) \gamma^5\sigma_3 \Xi  \eta(\vec y_2,t_0)G(\vec y_2,t_0; \vec x, t) ]

  //\vec p_1 + \vec p_2 = \vec p

  //\Xi is a projector 1/2(1+-sigma_2)

  //C(t) =   1/2 \sum_{\vec x} e^{i( \vec p\cdot \vec x ) } tr [ \gamma^5\sigma_3 R(\vec x,t; \vec p_1) \gamma^5\sigma_3 \Xi  S(\vec x,t; \vec p_2) ]

  //R(\vec x, t; \vec p_1) = \sum_{\vec y_1} e^{-i \vec p_1 \cdot \vec y_1} G(\vec x, t; \vec y_1, t_0)\eta(\vec y_1, t_0)

  //S(\vec x, t; \vec p_2) = \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2} \eta(\vec y_2, t_0)  G(\vec y_2, t_0; \vec x, t)
  //             = \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2} \eta(\vec y_2, t_0)\gamma^5 G^\dagger(\vec x, t; \vec y_2, t_0) \gamma^5
  //             = \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2} \eta(\vec y_2, t_0)\gamma^5 [ \gamma^5 C^{-1} \sigma_2 G^T(\vec x, t; \vec y_2, t_0) C\gamma^5 \sigma_2] \gamma^5
  //             = \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2} \eta(\vec y_2, t_0)C^{-1} \sigma_2 G^T(\vec x, t; \vec y_2, t_0) C \sigma_2

  //Iff \eta(\vec y_2, t_0) has diagonal spin, color, flavor structure

  //             = C^{-1} \sigma_2 [ \sum_{\vec y_2} e^{-i \vec p_2 \cdot \vec y_2} G(\vec x, t; \vec y_2, t_0)\eta(\vec y_2, t_0) ]^T C \sigma_2
  //             = C^{-1} \sigma_2 R(\vec x, t; \vec p_2)^T C \sigma_2

  //Should think about how this works with gauge fixing matrices in the sources themselves. Currently we just gauge fix the links themselves

  std::vector<RealD> momWallSourcePionCorrelator(const std::vector<int> &p1, const std::vector<int> &p2, int t0, const LatticeSCFmatrixD& R_p1, const LatticeSCFmatrixD& R_p2){
    std::cout << GridLogMessage << "Starting momentum wall source pion correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R_p1.Grid();
    std::vector<int> p = { p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2] };
    std::vector<double> pphys = getPhysicalMomentum(p);
    std::cout << GridLogMessage << "Total momentum in base units = " << p << " and physical units = " << pphys << std::endl;
    LatticeComplexD snk_phase_field = conjugate(phaseField(pphys, R_p1.Grid())); //exp(i \vec p \cdot \vec x)
  
    GparityFlavour sigma3(GparityFlavour::Algebra::SigmaZ);
    GparityFlavour sigma2(GparityFlavour::Algebra::SigmaY);
    Gamma gamma5(Gamma::Algebra::Gamma5);
    Gamma gamma2(Gamma::Algebra::GammaY);
    Gamma mgamma2(Gamma::Algebra::MinusGammaY);
    Gamma gamma4(Gamma::Algebra::MinusGammaT); //not GammaT?
  
    //C=-\gamma^2 \gamma^4
    //C^{-1} = -C
    LatticeSCFmatrixD S_p2 = gamma2 * ( gamma4 * ( sigma2 * transpose(R_p2) * mgamma2 ) * gamma4 )* sigma2;

    GparityFlavour Xi(getProjector(p2));
    std::cout << GridLogMessage << "Projector for p2=" << p2 << " : " << getProjector(p2) << std::endl;

    LatticeSCFmatrixD sqb = snk_phase_field * ( sigma3* ( gamma5 * R_p1 * gamma5 ) * sigma3 ) * Xi * S_p2;
  
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


  //pion propagator for regular action
  //C(t) =   \sum_{\vec x} tr [ \gamma^5 R(\vec x,t) \gamma^5  S(\vec x,t) ]

  //R(\vec x, t) = \sum_{\vec y_1} G(\vec x, t; \vec y_1, t_0)\eta(\vec y_1, t_0)

  //S(\vec x, t) = \sum_{\vec y_2} \eta(\vec y_2, t_0)  G(\vec y_2, t_0; \vec x, t)
  //             = \sum_{\vec y_2} \eta(\vec y_2, t_0)\gamma^5 G^\dagger(\vec x, t; \vec y_2, t_0) \gamma^5
  //Iff \eta(\vec y_2, t_0) has is real/hermitian and has diagonal spin structure
  //             = \gamma^5 R(\vec x,t)^\dagger \gamma^5

  std::vector<RealD> pionCorrelator(int t0, const LatticePropagatorD& R){
    std::cout << GridLogMessage << "Starting pion correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R.Grid();
  
    LatticePropagatorD sqb = R * adj(R) ;
  
    std::vector<SpinColourMatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    for(int t=0;t<Lt;t++)
      std::cout << GridLogMessage << "C(t)=" << trace(Ctm[t]) << std::endl;
    
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace(Ctm[t]) );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }
  //S is the backwards propagating quark
  std::vector<RealD> pionCorrelator(int t0, const LatticePropagatorD& R, const LatticePropagatorD &S){
    std::cout << GridLogMessage << "Starting pion correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R.Grid();
    Gamma gamma5(Gamma::Algebra::Gamma5);
    LatticePropagatorD sqb = gamma5 * (R * (gamma5 * S)) ;
  
    std::vector<SpinColourMatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    for(int t=0;t<Lt;t++)
      std::cout << GridLogMessage << "C(t)=" << trace(Ctm[t]) << std::endl;
    
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace(Ctm[t]) );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }


  std::vector<RealD> wallSinkPionCorrelator(int t0, const LatticePropagatorD& R){
    std::cout << GridLogMessage << "Starting wall sink pion correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R.Grid();
    
    std::vector<SpinColourMatrixD> Rsum;
    sliceSum(R, Rsum, 3);

    int Lt = Rsum.size();

    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace( Rsum[t]  * adj(Rsum[t]) ) );

    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }

  //pion to local axial-t propagator for regular action
  //C(t) =   \sum_{\vec x} tr [ \gamma^0\gamma^5 R(\vec x,t) \gamma^5  S(\vec x,t) ]

  //R(\vec x, t) = \sum_{\vec y_1} G(\vec x, t; \vec y_1, t_0)\eta(\vec y_1, t_0)

  //S(\vec x, t) = \sum_{\vec y_2} \eta(\vec y_2, t_0)  G(\vec y_2, t_0; \vec x, t)
  //             = \sum_{\vec y_2} \eta(\vec y_2, t_0)\gamma^5 G^\dagger(\vec x, t; \vec y_2, t_0) \gamma^5
  //Iff \eta(\vec y_2, t_0) has is real/hermitian and has diagonal spin structure
  //             = \gamma^5 R(\vec x,t)^\dagger \gamma^5

  std::vector<RealD> pionToLocalAxialTcorrelator(int t0, const LatticePropagatorD& R){
    std::cout << GridLogMessage << "Starting pion to local axial-t correlator with t0=" << t0 << std::endl;
    GridBase* UGrid = R.Grid();
    Gamma gamma4(Gamma::Algebra::GammaT);
    LatticePropagatorD sqb = gamma4 * (R * adj(R)) ;
  
    std::vector<SpinColourMatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    for(int t=0;t<Lt;t++)
      std::cout << GridLogMessage << "C(t)=" << trace(Ctm[t]) << std::endl;
    
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace(Ctm[t]) );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }
  //Requires 5D propagator solution
  template<typename ActionD> 
  std::vector<RealD> pionToConservedAxialTcorrelator(int t0, const LatticePropagatorD& R5D, ActionD &action){
    std::cout << GridLogMessage << "Starting pion to conserved axial-t correlator with t0=" << t0 << std::endl;
    assert(R5D.Grid() == action.FermionGrid());

    GridBase* UGrid = action.GaugeGrid();
    LatticePropagatorD fake_src(UGrid); fake_src = Zero(); //in principle you need the source to deal with the contact term, but we will never include data with tsep = 0

    LatticePropagatorD sqb(UGrid);
    LatticePropagatorD &Ru = const_cast<LatticePropagatorD &>(R5D); //grr!
    action.ContractConservedCurrent(Ru,Ru,sqb,fake_src,Current::Axial,3);
  
    std::vector<SpinColourMatrixD> Ctm;
    sliceSum(sqb, Ctm, 3);
  
    int Lt = Ctm.size();
    std::cout << GridLogMessage << "Computed correlator for " << Lt << " timeslices" << std::endl;

    for(int t=0;t<Lt;t++)
      std::cout << GridLogMessage << "C(t)=" << trace(Ctm[t]) << std::endl;
    
    std::vector<RealD> Ctx(Lt); //time coordinate is x[3]
    for(int t=0;t<Lt;t++)
      Ctx[t] = real( trace(Ctm[t]) );

    std::vector<RealD> Ct(Lt); //time coordinate is t = x[3] - y0[3]
    for(int tx=0;tx<Lt;tx++){
      int t = ( tx - t0 + Lt ) % Lt;
      Ct[t] = Ctx[tx];
    }
    return Ct;
  }
 
};
