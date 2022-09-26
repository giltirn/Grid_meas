#pragma once

#include "defines.h"

namespace GridMeas{
  using namespace Grid;

  //Convience package for the various common Grids
  struct Grids{
    //4d grids
    GridCartesian* UGrid;
    GridRedBlackCartesian* UrbGrid;
    
    //5d grids
    GridCartesian* FGrid;
    GridRedBlackCartesian* FrbGrid;

    int Ls;
    int subgrid_volume; //in the context of split grid, how many Grids make up a full Grid 
    
    Grids(): UGrid(nullptr), UrbGrid(nullptr), FGrid(nullptr), FrbGrid(nullptr), Ls(0), subgrid_volume(1){}
    ~Grids(){
      if(UGrid) delete UGrid;
      if(UrbGrid) delete UrbGrid;
      if(FGrid) delete FGrid;
      if(FrbGrid) delete FrbGrid;
    }
    
    void generate(const Coordinate &latt, const int _Ls, const int Nsimd){
      Ls = _Ls;
      subgrid_volume = 1;
      UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, GridDefaultSimd(Nd, Nsimd), GridDefaultMpi());
      UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
      FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
      FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
    }
    //Generate the Grids for split Grid
    //subgrid_proc_geom : the geometry of the MPI ranks of the subgrid
    void split(const Grids &r, const Coordinate &subgrid_proc_geom){
      Ls = r.Ls;
      UGrid = new GridCartesian(r.UGrid->_fdimensions,
				r.UGrid->_simd_layout,
				subgrid_proc_geom,
				*r.UGrid); 
      
      FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
      UrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
      FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
      subgrid_volume = 1;
      for(int i=0;i<Nd;i++)
	subgrid_volume *= r.UGrid->_processors[i]/subgrid_proc_geom[i];      
    }
    
  };
 
  static inline Grids makeDoublePrecGrids(const int Ls, const Coordinate &latt = GridDefaultLatt()){
    Grids out; out.generate(latt,Ls,vComplexD::Nsimd()); return out;
  }
  static inline Grids makeSinglePrecGrids(const int Ls, const Coordinate &latt = GridDefaultLatt()){
    Grids out; out.generate(latt,Ls,vComplexF::Nsimd()); return out;
  }
  //subgrid_proc_geom : the geometry of the MPI ranks of the subgrid
  //e.g. for an 4^4 original MPI layout generate 16 subgrids of 2^4 processors each with subgrid_proc_geom = (2,2,2,2)
  static inline Grids makeSplitGrids(const Grids &r, const Coordinate &subgrid_proc_geom){
    Grids out; out.split(r,subgrid_proc_geom); return out;
  }
    
    
      


};
