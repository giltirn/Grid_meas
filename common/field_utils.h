#pragma once

#include "defines.h"
#include "grids.h"
#include "field_utils.h"

namespace GridMeas{
  using namespace Grid;

  //From the matrices extract the column with the given flavor, spin and color index
  FermionFieldD extractColumn(const LatticeSCFmatrixD &from, const int fs, const int ss, const int cs){
    FermionFieldD into(from.Grid());
    {
      size_t Nsimd = FermionFieldD::vector_type::Nsimd();
      autoView( from_v, from, AcceleratorRead);
      autoView( into_v, into, AcceleratorWrite);
      accelerator_for( i, from_v.size(), Nsimd, {
	  auto site_matrix = from_v(i);
	  auto site_into = into_v(i);
	  for(int fr=0;fr<Ngp;fr++)
	    for(int sr=0;sr<Ns;sr++)
	      for(int cr=0;cr<Nc;cr++)
		site_into(fr)(sr)(cr) = site_matrix(fr,fs)(sr,ss)(cr,cs);
	  coalescedWrite(into_v[i], site_into);
	});
    }
    return into;
  }

  void insertColumn(LatticeSCFmatrixD &into, const FermionFieldD &from, const int fs, const int ss, const int cs){
    size_t Nsimd = FermionFieldD::vector_type::Nsimd();
    autoView( from_v, from, AcceleratorRead);
    autoView( into_v, into, AcceleratorWrite);
    accelerator_for( i, from_v.size(), Nsimd, {
	auto site_spinor = from_v(i);
	auto site_into = into_v(i);
	for(int fr=0;fr<Ngp;fr++)
	  for(int sr=0;sr<Ns;sr++)
	    for(int cr=0;cr<Nc;cr++)
	      site_into(fr,fs)(sr,ss)(cr,cs) = site_spinor(fr)(sr)(cr);
	coalescedWrite(into_v[i], site_into);
      });
  }

  //From the matrices extract the column with the given spin and color index
  FermionField1fD extractColumn(const LatticePropagatorD &from, const int ss, const int cs){
    FermionField1fD into(from.Grid());
    {
      size_t Nsimd = FermionField1fD::vector_type::Nsimd();
      autoView( from_v, from, AcceleratorRead);
      autoView( into_v, into, AcceleratorWrite);
      accelerator_for( i, from_v.size(), Nsimd, {
	  auto site_matrix = from_v(i);
	  auto site_into = into_v(i);
	    for(int sr=0;sr<Ns;sr++)
	      for(int cr=0;cr<Nc;cr++)
		site_into()(sr)(cr) = site_matrix()(sr,ss)(cr,cs);
	  coalescedWrite(into_v[i], site_into);
	});
    }
    return into;
  }

  void insertColumn(LatticePropagatorD &into, const FermionField1fD &from, const int ss, const int cs){
    size_t Nsimd = FermionField1fD::vector_type::Nsimd();
    autoView( from_v, from, AcceleratorRead);
    autoView( into_v, into, AcceleratorWrite);
    accelerator_for( i, from_v.size(), Nsimd, {
	auto site_spinor = from_v(i);
	auto site_into = into_v(i);
	  for(int sr=0;sr<Ns;sr++)
	    for(int cr=0;cr<Nc;cr++)
	      site_into()(sr,ss)(cr,cs) = site_spinor()(sr)(cr);
	coalescedWrite(into_v[i], site_into);
      });
  }


  template<typename FieldType>
  struct columnOps{};

  template<>
  struct columnOps<FermionFieldD>{
    typedef LatticeSCFmatrixD propagatorType;
    enum {Ncols=Ns*Nc*Ngp};
    inline static void unlex(int &f, int &s, int &c, const int col){
      int rem = col;
      c = rem % Nc; rem /= Nc;
      s = rem % Ns; rem /= Ns;
      f = rem;
    }
    inline static FermionFieldD extractColumn(const LatticeSCFmatrixD &from, const int col){
      int f,s,c; unlex(f,s,c,col); return GridMeas::extractColumn(from,f,s,c);
    }
    inline static void insertColumn(LatticeSCFmatrixD &into, const FermionFieldD &from, const int col){
      int f,s,c; unlex(f,s,c,col); return GridMeas::insertColumn(into,from,f,s,c);
    }
  };
  template<>
  struct columnOps<FermionField1fD>{
    typedef LatticePropagatorD propagatorType;
    enum {Ncols=Ns*Nc};
    inline static void unlex(int &s, int &c, const int col){
      int rem = col;
      c = rem % Nc; rem /= Nc;
      s = rem;
    }
    inline static FermionField1fD extractColumn(const LatticePropagatorD &from, const int col){
      int s,c; unlex(s,c,col); return GridMeas::extractColumn(from,s,c);
    }
    inline static void insertColumn(LatticePropagatorD &into, const FermionField1fD &from, const int col){
      int s,c; unlex(s,c,col); return GridMeas::insertColumn(into,from,s,c);
    }
  };

  


  //The midpoint propagator
  //\phi(x) = + P_L \psi(x,Ls/2) + P_R \psi(x,Ls/2-1) 
  template<typename Impl>
  typename Impl::FermionField extractMidProp(const typename Impl::FermionField &sol5d, CayleyFermion5D<Impl> &action_d){
    typedef typename Impl::FermionField Ferm;
    GridBase* UGrid = action_d.GaugeGrid();
    GridBase* FGrid = action_d.FermionGrid();
    
    Gamma G5(Gamma::Algebra::Gamma5);
    int Ls = action_d.Ls;
    Ferm p_plus (UGrid);			  
    Ferm p_minus(UGrid);
    Ferm p(UGrid);

    ExtractSlice(p_plus , sol5d, Ls/2-1 , 0);
    ExtractSlice(p_minus, sol5d, Ls/2   , 0);
    p_plus = p_plus + G5*p_plus;
    p_minus= p_minus - G5*p_minus;
    p=0.5*(p_plus+p_minus);
    return p;
  }



};
