#pragma once

#include "defines.h"
#include "grids.h"
#include "utils.h"

namespace GridMeas{
  using namespace Grid;

  //Template classes for picking the appropriate guesser
  template<typename FieldTypeD, typename EvecFieldType, int prec>
  struct _get_guesser{};

  //Create a guesser for a specific evec field type
  template<typename FieldTypeD, typename EvecFieldType>
  inline LinearFunction<FieldTypeD>* getGuesser(std::vector<Real> const& evals, std::vector<EvecFieldType> const &evecs){
    return _get_guesser<FieldTypeD, EvecFieldType, getPrecision<EvecFieldType>::value>::get(evals, evecs);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for double prec evecs
  template<typename FieldTypeD>
  struct _get_guesser<FieldTypeD, FieldTypeD, 2>{    
    inline static LinearFunction<FieldTypeD>* get(std::vector<Real> const& evals, std::vector<FieldTypeD> const &evecs){
      std::cout << GridLogMessage << "Using double precision guesser" << std::endl;
      return new DeflatedGuesser<FieldTypeD>(evecs, evals);
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for single prec evecs
  template<class FieldF, class FieldD>
  class MixedPrecDeflatedGuesser: public LinearFunction<FieldD> {
  private:
    const std::vector<FieldF> &evec;
    const std::vector<RealD> &eval;
    
  public:
    
    MixedPrecDeflatedGuesser(const std::vector<FieldF> & _evec,const std::vector<RealD> & _eval) : evec(_evec), eval(_eval) {};
    
    virtual void operator()(const FieldD &src,FieldD &guess) {
      assert(evec.size() > 0);
      assert(evec.size()==eval.size());
      auto N = evec.size();

#if 1
      FieldF src_f(evec[0].Grid()), guess_f(evec[0].Grid());
      precisionChange(src_f, src);
      guess_f = Zero();
      
      for (int i=0;i<N;i++) {
	const FieldF& tmp = evec[i];
	axpy(guess_f,TensorRemove(innerProduct(tmp,src_f)) / eval[i],tmp,guess_f);
      }
      precisionChange(guess, guess_f);
      guess.Checkerboard() = src.Checkerboard();
#else //note the below is more precise but costs a precision change for every evec
      FieldD tmp(src.Grid());
      guess = Zero();
      for (int i=0;i<N;i++) {
	precisionChange(tmp, evec[i]);
	axpy(guess, TensorRemove(innerProduct(tmp,src)) / eval[i],tmp,guess);
      }
      guess.Checkerboard() = src.Checkerboard();
#endif
    }

  };


  template<typename FieldTypeD, typename FieldTypeF>
  struct _get_guesser<FieldTypeD, FieldTypeF, 1>{
    inline static LinearFunction<FieldTypeD>* get(std::vector<Real> const& evals, std::vector<FieldTypeF> const &evecs){
      std::cout << GridLogMessage << "Using mixed precision guesser" << std::endl;
      return new MixedPrecDeflatedGuesser<FieldTypeF,FieldTypeD>(evecs, evals);
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for double prec X-conjugate evecs
  template<class Field>
  class XconjDeflatedGuesser: public LinearFunction<Field> {
  private:
    const std::vector<Field> &evec;
    const std::vector<RealD> &eval;
    const unsigned int       N;

  public:
    using LinearFunction<Field>::operator();

    XconjDeflatedGuesser(const std::vector<Field> & _evec,const std::vector<RealD> & _eval)
      : XconjDeflatedGuesser(_evec, _eval, _evec.size())
    {}

    XconjDeflatedGuesser(const std::vector<Field> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
      : evec(_evec), eval(_eval), N(_N)
    {
      assert(evec.size()==eval.size());
      assert(N <= evec.size());
    } 

    virtual void operator()(const Field &src,Field &guess) {
      guess = Zero();
      for (int i=0;i<N;i++) {
	const Field& tmp = evec[i];
	axpy(guess, 2.0*real(TensorRemove(innerProduct(tmp,src))) / eval[i],tmp,guess);
      }
      guess.Checkerboard() = src.Checkerboard();
    }
  };

  template<>
  struct _get_guesser<FermionField1fD, FermionField1fD, 2>{    
    inline static LinearFunction<FermionField1fD>* get(std::vector<Real> const& evals, std::vector<FermionField1fD> const &evecs){
      std::cout << GridLogMessage << "Using double precision Xconj guesser" << std::endl;
      return new XconjDeflatedGuesser<FermionField1fD>(evecs, evals);
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for single prec X-conjugate evecs
  template<class FieldF, class FieldD>
  class XconjMixedPrecDeflatedGuesser: public LinearFunction<FieldD> {
  private:
    const std::vector<FieldF> &evec;
    const std::vector<RealD> &eval;
    
  public:
    
    XconjMixedPrecDeflatedGuesser(const std::vector<FieldF> & _evec,const std::vector<RealD> & _eval) : evec(_evec), eval(_eval) {};
    
    virtual void operator()(const FieldD &src,FieldD &guess) {
      assert(evec.size() > 0);
      assert(evec.size()==eval.size());
      auto N = evec.size();
      FieldF src_f(evec[0].Grid()), guess_f(evec[0].Grid());
      precisionChange(src_f, src);
      guess_f = Zero();
      
      for (int i=0;i<N;i++) {
	const FieldF& tmp = evec[i];
	axpy(guess_f,2.0*real(TensorRemove(innerProduct(tmp,src_f))) / eval[i],tmp,guess_f);
      }
      precisionChange(guess, guess_f);
      guess.Checkerboard() = src.Checkerboard();
    }
  };

  template<>
  struct _get_guesser<FermionField1fD, FermionField1fF, 1>{
    inline static LinearFunction<FermionField1fD>* get(std::vector<Real> const& evals, std::vector<FermionField1fF> const &evecs){
      std::cout << GridLogMessage << "Using mixed precision Xconj guesser" << std::endl;
      return new XconjMixedPrecDeflatedGuesser<FermionField1fF,FermionField1fD>(evecs, evals);
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for double prec X-conjugate evecs with conversion to 2f fields, suitable for using X-conj evecs with regular GP Dirac op
  class XconjDeflatedGuesser2fConvert: public LinearFunction<FermionFieldD> {
  private:
    const std::vector<FermionField1fD> &evec;
    const std::vector<RealD> &eval;
    const unsigned int       N;
  public:
    using LinearFunction<FermionFieldD>::operator();

    XconjDeflatedGuesser2fConvert(const std::vector<FermionField1fD> & _evec,const std::vector<RealD> & _eval)
      : XconjDeflatedGuesser2fConvert(_evec, _eval, _evec.size())
    {}

    XconjDeflatedGuesser2fConvert(const std::vector<FermionField1fD> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
      : evec(_evec), eval(_eval), N(_N)
    {
      assert(evec.size()==eval.size());
      assert(N <= evec.size());
    } 

    virtual void operator()(const FermionFieldD &src,FermionFieldD &guess) {
      guess = Zero();
      FermionFieldD tmp2f(src.Grid());
      FermionField1fD tmp1f(src.Grid());

      for (int i=0;i<N;i++) {
	PokeIndex<GparityFlavourIndex>(tmp2f,evec[i],0);
	tmp1f = -(Xmatrix()*conjugate(evec[i]));
	PokeIndex<GparityFlavourIndex>(tmp2f,tmp1f,1);
	axpy(guess, real(TensorRemove(innerProduct(tmp2f,src))) / eval[i],tmp2f,guess);
      }
      guess.Checkerboard() = src.Checkerboard();
    }
  };
  template<>
  struct _get_guesser<FermionFieldD, FermionField1fD, 2>{    
    inline static LinearFunction<FermionFieldD>* get(std::vector<Real> const& evals, std::vector<FermionField1fD> const &evecs){
      std::cout << GridLogMessage << "Using double precision Xconj guesser with conversion to 2f" << std::endl;
      return new XconjDeflatedGuesser2fConvert(evecs, evals);
    }
  };


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Guesser for single prec X-conjugate evecs with conversion to 2f fields, suitable for using X-conj evecs with regular GP Dirac op
  class XconjMixedPrecDeflatedGuesser2fConvert: public LinearFunction<FermionFieldD> {
  private:
    const std::vector<FermionField1fF> &evec;
    const std::vector<RealD> &eval;
    const unsigned int       N;
  public:
    using LinearFunction<FermionFieldD>::operator();

    XconjMixedPrecDeflatedGuesser2fConvert(const std::vector<FermionField1fF> & _evec,const std::vector<RealD> & _eval)
      : XconjMixedPrecDeflatedGuesser2fConvert(_evec, _eval, _evec.size())
    {}

    XconjMixedPrecDeflatedGuesser2fConvert(const std::vector<FermionField1fF> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
      : evec(_evec), eval(_eval), N(_N)
    {
      assert(evec.size()==eval.size());
      assert(N <= evec.size());
    } 

    virtual void operator()(const FermionFieldD &src,FermionFieldD &guess) {
      FermionFieldF src_f(evec[0].Grid()), guess_f(evec[0].Grid());
      precisionChange(src_f, src);
      guess_f = Zero();
      FermionFieldF tmp2f(src_f.Grid());
      FermionField1fF tmp1f(src_f.Grid());

      for (int i=0;i<N;i++) {
	PokeIndex<GparityFlavourIndex>(tmp2f,evec[i],0);
	tmp1f = -(Xmatrix()*conjugate(evec[i]));
	PokeIndex<GparityFlavourIndex>(tmp2f,tmp1f,1);
	axpy(guess_f, real(TensorRemove(innerProduct(tmp2f,src_f))) / eval[i],tmp2f,guess_f);
      }

      precisionChange(guess, guess_f);
      guess.Checkerboard() = src.Checkerboard();
    }
  };
  template<>
  struct _get_guesser<FermionFieldD, FermionField1fF, 1>{    
    inline static LinearFunction<FermionFieldD>* get(std::vector<Real> const& evals, std::vector<FermionField1fF> const &evecs){
      std::cout << GridLogMessage << "Using single precision Xconj guesser with conversion to 2f" << std::endl;
      return new XconjMixedPrecDeflatedGuesser2fConvert(evecs, evals);
    }
  };


};
