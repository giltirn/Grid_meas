#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  inline void unlex(Coordinate &c, size_t p, const Coordinate &latt, int dims){
    for(int i=0;i<dims;i++){
      c[i] = p % latt[i];  p/=latt[i];
    }
  }
  
  bool operator==(const Coordinate &l, const Coordinate &r){
    if(l.size() != r.size()) return false;
    for(int i=0;i<l.size();i++)
      if(l[i] != r[i])
	return false;
    return true;
  }
  bool operator!=(const Coordinate &l, const Coordinate &r){ return !(l == r); }


  //Fill the first 'dims' elements of the coordinate 'c' with random values
  template<typename RNG>
  void randomCoordinate(Coordinate &c, const Coordinate &latt, int dims, RNG &std_rng, GridBase* any_grid){
    size_t vol = 1;
    for(int i=0;i<dims;i++) vol *= latt[i];

    std::uniform_int_distribution<size_t> distribution(0,vol-1);
    c.resize(Nd);
    size_t p = distribution(std_rng);
    any_grid->Broadcast(0, (void*)&p, sizeof(size_t));
    unlex(c, p, latt, dims);
  }



  //Get an instance of GetType from among 2 objects (type must match A or B)
  template<typename GetType, typename A, typename B>
  struct _getInstance2;

  template<typename A, typename B>
  struct _getInstance2<A,A,B>{
    static inline A & doit(A &a, B& b){ return a; }
  };
  template<typename A, typename B>
  struct _getInstance2<B,A,B>{
    static inline B & doit(A &a, B& b){ return b; }
  };

  template<typename GetType, typename A, typename B>
  inline GetType & getInstance(A &a, B& b){ return _getInstance2<GetType,A,B>::doit(a,b); }


  //Get an instance of GetType from among 4 objects (type must match A,B,C or D)
  template<typename GetType, typename A, typename B, typename C, typename D>
  struct _getInstance4;

  template<typename A, typename B, typename C, typename D>
  struct _getInstance4<A,A,B,C,D>{
    static inline A & doit(A &a, B& b, C& c, D& d){ return a; }
  };
  template<typename A, typename B, typename C, typename D>
  struct _getInstance4<B,A,B,C,D>{
    static inline B & doit(A &a, B& b, C& c, D& d){ return b; }
  };
  template<typename A, typename B, typename C, typename D>
  struct _getInstance4<C,A,B,C,D>{
    static inline C & doit(A &a, B& b, C& c, D& d){ return c; }
  };
  template<typename A, typename B, typename C, typename D>
  struct _getInstance4<D,A,B,C,D>{
    static inline D & doit(A &a, B& b, C& c, D& d){ return d; }
  };

  template<typename GetType, typename A, typename B, typename C, typename D>
  inline GetType & getInstance(A &a, B& b, C& c, D& d){ return _getInstance4<GetType,A,B,C,D>::doit(a,b,c,d); }

  inline Gamma Xmatrix(){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    return X;
  }

  inline void addResult(std::vector<RealD> &into, const std::vector<RealD> &from, const RealD nsrc){
    for(int t=0;t<into.size();t++) //average over sources
      into[t] += from[t] / nsrc;
  }

  inline void addResult(std::vector<RealD> &into_avg,  std::vector<std::vector<RealD> > &into_sep, const std::vector<RealD> &from, const RealD nsrc){
    for(int t=0;t<into_avg.size();t++) //average over sources
      into_avg[t] += from[t] / nsrc;

    into_sep.push_back(from);
  }


}
