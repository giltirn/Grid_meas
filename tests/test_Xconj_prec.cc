#include<common/propagator_invert_field_Xconj.h>
#include<common/propagator_invert_field.h>
#include<common/utils.h>
#include<common/action.h>
#include<common/sources.h>
#include<common/field_utils.h>

using namespace Grid;
using namespace GridMeas;

template<typename Op>
void applyM(FermionField1fD &out, const FermionField1fD &in, const int f1, const int f2, CayleyFermion5D<GparityWilsonImplD> &action, const Op &op){
  FermionFieldD in2f(in.Grid()), out2f(in.Grid());
  in2f.Checkerboard() = in.Checkerboard();
  in2f = Zero();
  PokeIndex<GparityFlavourIndex>(in2f, in, f2);
  op(out2f, in2f, action);
  out = PeekIndex<GparityFlavourIndex>(out2f,f1);
  out.Checkerboard() = out2f.Checkerboard();
}

//Note inv only applies if cbin=cbout
template<typename Action>
bool applyOp(typename Action::FermionField &out, typename Action::FermionField &in, int cbout, int cbin, int dag, int inv, Action &action){
  in.Checkerboard() = cbin;
  out.Checkerboard() = cbout;

  if(cbin == cbout){
    if(dag){
      if(inv) action.MooeeInvDag(in,out);
      else    action.MooeeDag(in,out);
    }else{
      if(inv) action.MooeeInv(in,out);
      else    action.Mooee(in,out);
    }
  }else{
    if(inv) return false;

    if(dag) action.MeooeDag(in,out);
    else    action.Meooe(in,out);
  }
  return true;
}


int main(int argc, char** argv){
  Grid_init(&argc, &argv);

  int Ls = 8;
  Coordinate latt = GridDefaultLatt();
  int Lt = latt[3];
  size_t V4d = latt[0]*latt[1]*latt[2]*latt[3];
  size_t V3d = latt[0]*latt[1]*latt[2];

  Grids GridsD = makeDoublePrecGrids(Ls, latt);
  Grids GridsF = makeSinglePrecGrids(Ls, latt);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG pRNG(GridsD.UGrid); //4D!
  pRNG.SeedFixedIntegers(seeds4);

  GridParallelRNG pRNG5(GridsD.FGrid);  
  pRNG5.SeedFixedIntegers(seeds4);

  GridSerialRNG sRNG;  
  sRNG.SeedFixedIntegers(seeds4); 

  LatticeGaugeFieldD Ud(GridsD.UGrid);
  LatticeGaugeFieldF Uf(GridsF.UGrid);
  SU<Nc>::HotConfiguration(pRNG, Ud);
  precisionChange(Uf,Ud);
  
  std::vector<int> dirs4({1,1,1,0});
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  Params.twists = Coordinate(std::vector<int>({1,1,1,1})); //APBC in t direction

  Actions actions(ActionType::Mobius, Params, 0.01, 2.0, Ud, GridsD, Uf, GridsF); 

  //Generate an X-conjugate random field
  FermionField1fD src_1f(GridsD.FrbGrid);
  gaussian(pRNG,src_1f);
  FermionFieldD src_2f(GridsD.FrbGrid);
  get2fXconjVector(src_2f, src_1f);

  std::cout << "Test full (uncheckerboarded) operation action on X-conjugate vectors from G-parity action components" << std::endl;
  { 
    FermionField1fD src_1f_full(GridsD.FGrid);
    gaussian(pRNG,src_1f_full);

    FermionField1fD result_1f(GridsD.FGrid);
    actions.xconj_action_d->M(src_1f_full, result_1f);

    FermionField1fD M11a(GridsD.FGrid);
    applyM(M11a, src_1f_full, 0,0, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.M(i,o); });

    FermionField1fD tmp2(GridsD.FGrid), mM12Xastar(GridsD.FGrid);
    tmp2 = -(Xmatrix()*conjugate(src_1f_full));
    applyM(mM12Xastar, tmp2, 0,1, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.M(i,o); });

    FermionField1fD result_2f = M11a + mM12Xastar;

    tmp2 = result_2f - result_1f;
    std::cout << "M_X a = M_11 a - M_12 Xa* :   lhs:" << norm2(result_1f) << " rhs:" << norm2(result_2f) << " diff (expect 0):" << norm2(tmp2) << std::endl;
    assert(norm2(tmp2) < 1e-10);
  }


  std::cout << "Test half-checkerboard operation action on X-conjugate vectors from G-parity action components" << std::endl;
  { 
    FermionField1fD result_1f(GridsD.FrbGrid);
    actions.xconj_action_d->Meooe(src_1f, result_1f);

    FermionField1fD M11a(GridsD.FrbGrid);
    applyM(M11a, src_1f, 0,0, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.Meooe(i,o); });

    FermionField1fD tmp2(GridsD.FrbGrid), mM12Xastar(GridsD.FrbGrid);
    tmp2 = -(Xmatrix()*conjugate(src_1f));
    applyM(mM12Xastar, tmp2, 0,1, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.Meooe(i,o); });

    FermionField1fD result_2f = M11a + mM12Xastar;

    tmp2 = result_2f - result_1f;
    std::cout << "For eo/oe    M_X a = M_11 a - M_12 Xa* :   lhs:" << norm2(result_1f) << " rhs:" << norm2(result_2f) << " diff (expect 0):" << norm2(tmp2) << std::endl;
    assert(norm2(tmp2) < 1e-10);
  }
  { 
    FermionField1fD result_1f(GridsD.FrbGrid);
    actions.xconj_action_d->Mooee(src_1f, result_1f);

    FermionField1fD M11a(GridsD.FrbGrid);
    applyM(M11a, src_1f, 0,0, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.Mooee(i,o); });

    FermionField1fD tmp2(GridsD.FrbGrid), mM12Xastar(GridsD.FrbGrid);
    tmp2 = -(Xmatrix()*conjugate(src_1f));
    applyM(mM12Xastar, tmp2, 0,1, *actions.action_d, [](FermionFieldD &o, const FermionFieldD &i, CayleyFermion5D<GparityWilsonImplD> &action){ action.Mooee(i,o); });

    FermionField1fD result_2f = M11a + mM12Xastar;

    tmp2 = result_2f - result_1f;
    std::cout << "For ee/oo    M_X a = M_11 a - M_12 Xa* :   lhs:" << norm2(result_1f) << " rhs:" << norm2(result_2f) << " diff (expect 0):" << norm2(tmp2) << std::endl;
    assert(norm2(tmp2) < 1e-10);
  }

  std::cout << "Testing half-checkerboard operation action on X-conjugate vectors" << std::endl;

  for(int cbout=0;cbout<=1;cbout++){
    for(int cbin=0;cbin<=1;cbin++){
      for(int dag=0;dag<=1;dag++){
	for(int inv=0;inv<=1;inv++){
	  FermionFieldD out_2f(GridsD.FrbGrid);
	  bool r2f = applyOp(out_2f, src_2f, cbout, cbin, dag, inv, *actions.action_d);
	  FermionField1fD out_1f(GridsD.FrbGrid);
	  bool r1f = applyOp(out_1f, src_1f, cbout, cbin, dag, inv, *actions.xconj_action_d);	  

	  if(r1f && r2f){
	    FermionFieldD out_1f_conv2f(GridsD.FrbGrid);
	    get2fXconjVector(out_1f_conv2f, out_1f);
	    
	    FermionFieldD tmp = out_1f_conv2f - out_2f;
	    std::cout << "cbout:" << cbout << " cbin:" << cbin << " dag:" << dag << " inv:" << inv << " result GP:" << norm2(out_2f) << " result Xconj:" << norm2(out_1f_conv2f) << " diff (expect 0): " << norm2(tmp) << std::endl;
	  }
	}
      }
    }
  }
  
  SchurOperatorBase<FermionFieldD>* schurop_2f[2] = { new SchurDiagMooeeOperator<CayleyFermion5D<GparityWilsonImplD>,FermionFieldD>(*actions.action_d),
						      new SchurDiagOneOperator<CayleyFermion5D<GparityWilsonImplD>,FermionFieldD>(*actions.action_d) };

  SchurOperatorBase<FermionField1fD>* schurop_1f[2] = { new SchurDiagMooeeOperator<CayleyFermion5D<XconjugateWilsonImplD>,FermionField1fD>(*actions.xconj_action_d),
							new SchurDiagOneOperator<CayleyFermion5D<XconjugateWilsonImplD>,FermionField1fD>(*actions.xconj_action_d) };
  std::cout << "Testing full preconditioned operation action on X-conjugate vectors" << std::endl;
  
  for(int op=0;op<=1;op++){
    for(int dag=0;dag<=1;dag++){
      FermionFieldD out_2f(GridsD.FrbGrid);
      if(dag) schurop_2f[op]->MpcDag(src_2f,out_2f);
      else schurop_2f[op]->Mpc(src_2f,out_2f);
      
      FermionField1fD out_1f(GridsD.FrbGrid);
      if(dag) schurop_1f[op]->MpcDag(src_1f,out_1f);
      else schurop_1f[op]->Mpc(src_1f,out_1f);

      FermionFieldD out_1f_conv2f(GridsD.FrbGrid);
      get2fXconjVector(out_1f_conv2f, out_1f);
      
      FermionFieldD tmp = out_1f_conv2f - out_2f;

      std::cout << (op == 0 ? "SchurOriginal" : "SchurDiagOne") << " dag:" << dag << " result GP:" << norm2(out_2f) << " result Xconj:" << norm2(out_1f_conv2f) << " diff (expect 0): " << norm2(tmp) << std::endl;
    
    }
  }

  std::cout << "Testing full preconditioned squared Hermitian operation action on X-conjugate vectors" << std::endl;
  
  for(int op=0;op<=1;op++){
    FermionFieldD out_2f(GridsD.FrbGrid);
    schurop_2f[op]->HermOp(src_2f,out_2f);
    
    FermionField1fD out_1f(GridsD.FrbGrid);
    schurop_1f[op]->HermOp(src_1f,out_1f);
    
    FermionFieldD out_1f_conv2f(GridsD.FrbGrid);
    get2fXconjVector(out_1f_conv2f, out_1f);
    
    FermionFieldD tmp = out_1f_conv2f - out_2f;
    
    std::cout << (op == 0 ? "SchurOriginal" : "SchurDiagOne") << " result GP:" << norm2(out_2f) << " result Xconj:" << norm2(out_1f_conv2f) << " diff (expect 0): " << norm2(tmp) << std::endl;
    
  }


  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
