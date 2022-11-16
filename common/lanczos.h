#pragma once

#include <Grid/Grid.h>
#include "array_io.h"
#include "field_array_io.h"

namespace GridMeas{
  using namespace Grid;

  struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
				    double, alpha,
				    double, beta,
				    double, mu,
				    int, ord,
				    int, n_stop,
				    int, n_want,
				    int, n_use,
				    double, tolerance);

    LanczosParameters() {
      alpha = 35;
      beta = 5;
      mu = 0;
      ord = 100;
      n_stop = 10;
      n_want = 10;
      n_use = 15;
      tolerance = 1e-6;
    }
  };

  template<typename ActionParams, typename Field>
  struct getLanczosInnerProduct{
    static void doit(innerProductImplementation<Field> *& inner, ImplicitlyRestartedLanczosTester<Field> *&tester, LinearFunction<Field> &HermOp){
      inner = new innerProductImplementation<Field>;
      tester = new ImplicitlyRestartedLanczosHermOpTester<Field>(HermOp, *inner);
    }
  };
  template<typename Field>
  struct getLanczosInnerProduct<XconjWilsonImplParams, Field>{
    static void doit(innerProductImplementation<Field> *& inner, ImplicitlyRestartedLanczosTester<Field> *&tester, LinearFunction<Field> &HermOp){
      inner = new innerProductImplementationXconjugate<Field>;
      tester = new ImplicitlyRestartedLanczosHermOpTester<Field>(HermOp, *inner);
    }
  };
   

  template<typename FermionAction, typename FermionField, typename LatticeGaugeField>
  void computeEigenvalues(std::vector<RealD> &eval,
			  std::vector<FermionField> &evec,			
			  const LanczosParameters &params,
			  GridCartesian* Grid, GridRedBlackCartesian* rbGrid, const LatticeGaugeField &latt,  //expect lattice to have been initialized to something
			  FermionAction &action, GridParallelRNG &rng){
    if(params.n_stop == 0 || params.n_want == 0){
      eval.resize(0); evec.resize(0,rbGrid);
      std::cout << "Skipping Lanczos as #evals required is zero" << std::endl;	
      return;
    }	
    
    FermionField gauss_o(rbGrid);
    FermionField gauss(Grid);
    gaussian(rng, gauss);
    pickCheckerboard(Odd, gauss_o, gauss);

    action.ImportGauge(latt);

    SchurDiagMooeeOperator<FermionAction, FermionField> hermop(action);
    PlainHermOp<FermionField> hermop_wrap(hermop);
    //ChebyshevLanczos<FermionField> Cheb(params.alpha, params.beta, params.mu, params.ord);
    assert(params.mu == 0.0);

    Chebyshev<FermionField> Cheb(params.beta*params.beta, params.alpha*params.alpha, params.ord+1);
    FunctionHermOp<FermionField> Cheb_wrap(Cheb, hermop);

    innerProductImplementation<FermionField> *inner;
    ImplicitlyRestartedLanczosTester<FermionField> *tester;
    getLanczosInnerProduct<typename FermionAction::ImplParams, FermionField>::doit(inner,tester,hermop_wrap);

    std::cout << "IRL: alpha=" << params.alpha << " beta=" << params.beta << " mu=" << params.mu << " ord=" << params.ord << std::endl;
    ImplicitlyRestartedLanczos<FermionField> IRL(Cheb_wrap, hermop_wrap, *tester, *inner, params.n_stop, params.n_want, params.n_use, params.tolerance, 10000);

    eval.resize(params.n_use);
    evec.clear();
    evec.resize(params.n_use, rbGrid);

    int Nconv;
    IRL.calc(eval, evec, gauss_o, Nconv);

    delete inner;
    delete tester;

    std::cout << "Converged " << eval.size() << " eigenvalues:" << std::endl;
    for(int i=0;i<eval.size();i++){
      std::cout << i << " " << eval[i] << std::endl;
    }
  }

  template<typename FermionField>
  void saveEigenvalues(const std::vector<RealD> &eval,
		       const std::vector<FermionField> &evec,
		       const std::string &eval_file_stub,
		       const std::string &evec_file_stub,
		       const int cfg){
    std::stringstream eval_file; eval_file << eval_file_stub << "." << cfg;
    std::stringstream evec_file; evec_file << evec_file_stub << "." << cfg;
    binaryWriteArray(eval_file.str(), eval);
    writeFieldArray(evec_file.str(), evec);
  }

  template<typename FermionField>
  void readEigenvalues(std::vector<RealD> &eval,
		       std::vector<FermionField> &evec,
		       GridBase* FrbGrid,
		       const std::string &eval_file_stub,
		       const std::string &evec_file_stub,
		       const int cfg){
    std::stringstream eval_file; eval_file << eval_file_stub << "." << cfg;
    std::stringstream evec_file; evec_file << evec_file_stub << "." << cfg;
    binaryReadArray(eval, eval_file.str());
    std::cout << GridLogMessage << "Read " << eval.size() << " eigenvalues" << std::endl;
    evec.resize(eval.size(), FrbGrid);
    readFieldArray(evec, evec_file.str());
    std::cout << GridLogMessage << "Read " << eval.size() << " eigenvectors" << std::endl;
  }

}
