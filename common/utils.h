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
}
