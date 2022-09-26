#pragma once

#include<Grid/Grid.h>

namespace GridMeas{
  using namespace Grid;

  typedef iMatrix<iMatrix<iMatrix<vComplexD, Nc>, Ns>, Ngp> vSCFmatrixD;
  typedef iMatrix<iMatrix<iMatrix<ComplexD, Nc>, Ns>, Ngp> SCFmatrixD;
  typedef Lattice<vSCFmatrixD> LatticeSCFmatrixD;

  typedef GparityWilsonImplD::FermionField FermionFieldD;
  typedef GparityWilsonImplD::SiteSpinor SiteSpinorD;
  
  typedef GparityWilsonImplF::FermionField FermionFieldF;
  typedef GparityWilsonImplF::SiteSpinor SiteSpinorF;

  typedef XconjugateWilsonImplD::FermionField FermionField1fD;
  typedef XconjugateWilsonImplD::SiteSpinor SiteSpinor1fD;
  
  typedef XconjugateWilsonImplF::FermionField FermionField1fF;
  typedef XconjugateWilsonImplF::SiteSpinor SiteSpinor1fF;

};
