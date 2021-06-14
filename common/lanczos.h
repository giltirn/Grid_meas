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



  template<typename FermionActionD, typename FermionFieldD>
  void computeEigenvalues(std::vector<RealD> &eval,
			  std::vector<FermionFieldD> &evec,			
			  const LanczosParameters &params,
			  GridCartesian* Grid, GridRedBlackCartesian* rbGrid, const LatticeGaugeFieldD &latt,  //expect lattice to have been initialized to something
			  FermionActionD &action, GridParallelRNG &rng){
  
    FermionFieldD gauss_o(rbGrid);
    FermionFieldD gauss(Grid);
    gaussian(rng, gauss);
    pickCheckerboard(Odd, gauss_o, gauss);

    action.ImportGauge(latt);

    SchurDiagMooeeOperator<FermionActionD, FermionFieldD> hermop(action);
    PlainHermOp<FermionFieldD> hermop_wrap(hermop);
    //ChebyshevLanczos<FermionFieldD> Cheb(params.alpha, params.beta, params.mu, params.ord);
    assert(params.mu == 0.0);

    Chebyshev<FermionFieldD> Cheb(params.beta*params.beta, params.alpha*params.alpha, params.ord+1);
    FunctionHermOp<FermionFieldD> Cheb_wrap(Cheb, hermop);

    std::cout << "IRL: alpha=" << params.alpha << " beta=" << params.beta << " mu=" << params.mu << " ord=" << params.ord << std::endl;
    ImplicitlyRestartedLanczos<FermionFieldD> IRL(Cheb_wrap, hermop_wrap, params.n_stop, params.n_want, params.n_use, params.tolerance, 10000);

    eval.resize(params.n_use);
    evec.clear();
    evec.resize(params.n_use, rbGrid);

    int Nconv;
    IRL.calc(eval, evec, gauss_o, Nconv);

    std::cout << "Eigenvalues:" << std::endl;
    for(int i=0;i<params.n_want;i++){
      std::cout << i << " " << eval[i] << std::endl;
    }
  }

  template<typename FermionFieldD>
  void saveEigenvalues(const std::vector<RealD> &eval,
		       const std::vector<FermionFieldD> &evec,
		       const std::string &eval_file_stub,
		       const std::string &evec_file_stub,
		       const int cfg){
    std::stringstream eval_file; eval_file << eval_file_stub << "." << cfg;
    std::stringstream evec_file; evec_file << evec_file_stub << "." << cfg;
    binaryWriteArray(eval_file.str(), eval);
    writeFieldArray(evec_file.str(), evec);
  }

  template<typename FermionFieldD>
  void readEigenvalues(std::vector<RealD> &eval,
		       std::vector<FermionFieldD> &evec,
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
