#pragma once

#include "defines.h"
#include "grids.h"

namespace GridMeas{
  using namespace Grid;

  GRID_SERIALIZABLE_ENUM(ActionType, undef, DWF, 1, Mobius, 2);

  CayleyFermion5D<GparityWilsonImplD>* createActionD(ActionType action,
						     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						     LatticeGaugeFieldD &Umu,
						     GridCartesian         &FiveDimGrid,
						     GridRedBlackCartesian &FiveDimRedBlackGrid,
						     GridCartesian         &FourDimGrid,
						     GridRedBlackCartesian &FourDimRedBlackGrid
						     ){
    double bpc = mobius_scale;
    double bmc = 1.0;
    double b = (bpc + bmc)/2.;
    double c = (bpc - bmc)/2.;
    
    switch(action){
    case ActionType::DWF:
      std::cout << GridLogMessage << "Creating double prec DWF action with m=" << mass << std::endl;
      return new GparityDomainWallFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, Params);
    case ActionType::Mobius:
      std::cout << GridLogMessage << "Creating double prec Mobius action with b+c=" << b+c << " b-c=" << b-c << " m=" << mass << std::endl;
      return new GparityMobiusFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, b, c, Params);
    };
    return nullptr;
  };

  CayleyFermion5D<GparityWilsonImplF>* createActionF(ActionType action,
						     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						     LatticeGaugeFieldF &Umu,
						     GridCartesian         &FiveDimGrid,
						     GridRedBlackCartesian &FiveDimRedBlackGrid,
						     GridCartesian         &FourDimGrid,
						     GridRedBlackCartesian &FourDimRedBlackGrid
						     ){
    double bpc = mobius_scale;
    double bmc = 1.0;
    double b = (bpc + bmc)/2.;
    double c = (bpc - bmc)/2.;
    
    switch(action){
    case ActionType::DWF:
      std::cout << GridLogMessage << "Creating single prec DWF action with m=" << mass << std::endl;	    
      return new GparityDomainWallFermionF(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, Params);
    case ActionType::Mobius:
      std::cout << GridLogMessage << "Creating single prec Mobius action with b+c=" << b+c << " b-c=" << b-c << " m=" << mass << std::endl;	    
      return new GparityMobiusFermionF(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, b, c, Params);
    };
    return nullptr;
  };


  CayleyFermion5D<GparityWilsonImplD>* createActionD(ActionType action,
						     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						     LatticeGaugeFieldD &Umu, Grids &Grids){
    return createActionD(action,Params,mass,mobius_scale,Umu,*Grids.FGrid,*Grids.FrbGrid,*Grids.UGrid,*Grids.UrbGrid);
  }
  CayleyFermion5D<GparityWilsonImplF>* createActionF(ActionType action,
						     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						     LatticeGaugeFieldF &Umu, Grids &Grids){
    return createActionF(action,Params,mass,mobius_scale,Umu,*Grids.FGrid,*Grids.FrbGrid,*Grids.UGrid,*Grids.UrbGrid);
  }



}
