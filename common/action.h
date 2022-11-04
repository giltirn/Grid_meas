#pragma once

#include "defines.h"
#include "grids.h"
#include "memory.h"

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


  XconjWilsonImplParams XconjParamsTransfer(const GparityWilsonImplParams &gpimpl){
    XconjWilsonImplParams out;
    out.dirichlet = gpimpl.dirichlet;
    out.twists = gpimpl.twists;
    out.boundary_phase = 1.0; //X-conj
    return out;
  }

  CayleyFermion5D<XconjugateWilsonImplD>* createXconjActionD(ActionType action,
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
    
    auto Xparams = XconjParamsTransfer(Params);

    switch(action){
    case ActionType::DWF:
      std::cout << GridLogMessage << "Creating double prec X-conjugate DWF action with m=" << mass << std::endl;
      return new XconjugateDomainWallFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, Xparams);
    case ActionType::Mobius:
      std::cout << GridLogMessage << "Creating double prec X-conjugate Mobius action with b+c=" << b+c << " b-c=" << b-c << " m=" << mass << std::endl;
      return new XconjugateMobiusFermionD(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, b, c, Xparams);
    };
    return nullptr;
  };

  CayleyFermion5D<XconjugateWilsonImplF>* createXconjActionF(ActionType action,
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
    
    auto Xparams = XconjParamsTransfer(Params);

    switch(action){
    case ActionType::DWF:
      std::cout << GridLogMessage << "Creating single prec X-conjugate DWF action with m=" << mass << std::endl;
      return new XconjugateDomainWallFermionF(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, Xparams);
    case ActionType::Mobius:
      std::cout << GridLogMessage << "Creating single prec X-conjugate Mobius action with b+c=" << b+c << " b-c=" << b-c << " m=" << mass << std::endl;
      return new XconjugateMobiusFermionF(Umu, FiveDimGrid, FiveDimRedBlackGrid, FourDimGrid, FourDimRedBlackGrid, mass, 1.8, b, c, Xparams);
    };
    return nullptr;
  };

  CayleyFermion5D<XconjugateWilsonImplD>* createXconjActionD(ActionType action,
							     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
							     LatticeGaugeFieldD &Umu, Grids &Grids){
    return createXconjActionD(action,Params,mass,mobius_scale,Umu,*Grids.FGrid,*Grids.FrbGrid,*Grids.UGrid,*Grids.UrbGrid);
  }
  CayleyFermion5D<XconjugateWilsonImplF>* createXconjActionF(ActionType action,
						     const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
						     LatticeGaugeFieldF &Umu, Grids &Grids){
    return createXconjActionF(action,Params,mass,mobius_scale,Umu,*Grids.FGrid,*Grids.FrbGrid,*Grids.UGrid,*Grids.UrbGrid);
  }

  struct Actions{
    CayleyFermion5D<GparityWilsonImplD>* action_d;
    CayleyFermion5D<GparityWilsonImplF>* action_f;
    CayleyFermion5D<XconjugateWilsonImplD>* xconj_action_d;
    CayleyFermion5D<XconjugateWilsonImplF>* xconj_action_f;

  
    Actions(): action_d(nullptr), action_f(nullptr), xconj_action_d(nullptr), xconj_action_f(nullptr){}

    Actions(ActionType action,
	    const GparityWilsonImplD::ImplParams &Params, double mass, double mobius_scale,
	    LatticeGaugeFieldD &U_d, Grids &GridsD,
	    LatticeGaugeFieldF &U_f, Grids &GridsF){
      printMem("Pre actionD creation");
      action_d = createActionD(action, Params, mass, mobius_scale, U_d, GridsD);
      printMem("Pre actionF creation");
      action_f = createActionF(action, Params, mass, mobius_scale, U_f, GridsF);
      printMem("Post actionD/F creation");

      printMem("Pre Xconj actionD creation");
      xconj_action_d = createXconjActionD(action, Params, mass, mobius_scale, U_d, GridsD);
      printMem("Pre Xconj actionF creation");
      xconj_action_f = createXconjActionF(action, Params, mass, mobius_scale, U_f, GridsF);
      printMem("Post Xconj actionD/F creation");
    }    
  
    template<typename ActionType>
    inline ActionType* getAction(){ return getInstance<ActionType*>(action_d, action_f, xconj_action_d, xconj_action_f); }

    ~Actions(){
      if(action_d) delete action_d;
      if(action_f) delete action_f;
      if(xconj_action_d) delete xconj_action_d;
      if(xconj_action_f) delete xconj_action_f;
    }

    void ImportGauge(LatticeGaugeFieldD &U_d,LatticeGaugeFieldF &U_f){
      if(action_d) action_d->ImportGauge(U_d);
      if(action_f) action_f->ImportGauge(U_f);
      if(xconj_action_d) xconj_action_d->ImportGauge(U_d);
      if(xconj_action_f) xconj_action_f->ImportGauge(U_f);
    }
		   
    
  };

  //Setup the gauge BCs and return the parameters to set up the fermion actions
  //Antiperiodic BCs in the time direction are assumed
  inline GparityWilsonImplD::ImplParams setupActionParams(const std::vector<Integer> &GparityDirs){
    assert(Nd == 4);
    std::vector<int> dirs4(4);
    for(int i=0;i<3;i++) dirs4[i] = GparityDirs[i];
    dirs4[3] = 0; //periodic gauge BC in time
    
    std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
    ConjugateGimplD::setDirections(dirs4); //gauge BC
    
    GparityWilsonImplD::ImplParams Params;
    for(int i=0;i<Nd-1;i++) Params.twists[i] = GparityDirs[i]; //G-parity directions
    Params.twists[Nd-1] = 1; //APBC in time direction
    std::cout << GridLogMessage << "Fermion BCs: " << Params.twists << std::endl;
    return Params;
  }

}
