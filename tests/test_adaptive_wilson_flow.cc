#include<common/wilson_flow.h>
#include<common/utils.h>
#include<common/grids.h>

using namespace Grid;
using namespace GridMeas;

//Linearly interpolate between two nearest times
RealD interpolate(const RealD t_int, const std::vector<std::pair<RealD,RealD> > &data){
  RealD tdiff1=1e32; int t1_idx=-1;
  RealD tdiff2=1e32; int t2_idx=-1;

  for(int i=0;i<data.size();i++){
    RealD diff = fabs(data[i].first-t_int);
    //std::cout << "targ " << t_int << " cur " << data[i].first << " diff " << diff << " best diff1 " << tdiff1 << " diff2 " << tdiff2 << std::endl;

    if(diff < tdiff1){ 
      if(tdiff1 < tdiff2){ //swap out tdiff2
	tdiff2 = tdiff1; t2_idx = t1_idx;
      }
      tdiff1 = diff; t1_idx = i; 
    }
    else if(diff < tdiff2){ tdiff2 = diff; t2_idx = i; }
  }
  assert(t1_idx != -1 && t2_idx != -1);
  
  RealD t2 = data[t2_idx].first,  v2 = data[t2_idx].second;
  RealD t1 = data[t1_idx].first,  v1 = data[t1_idx].second;
  
  //v = a + bt
  //v2-v1 = b(t2-t1)
  RealD b = (v2-v1)/(t2-t1);
  RealD a = v1 - b*t1;
  RealD vout = a + b*t_int;

  //std::cout << "Interpolate to " << t_int << " two closest points " << t1  << " " << t2 
  //<< " with values " << v1 << " "<< v2 << " : got " << vout <<  std::endl;
  return vout;
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

  WilsonFlowIO io_base;
  io_base.do_energy_density_clover = true;

  WilsonFlowIO io_orig = io_base;
  WilsonFlowMeasGeneral(io_orig, 300, 0.01, Ud);

  WilsonFlowIO io_adaptive = io_base;
  WilsonFlowAdaptiveMeasGeneral(io_adaptive, 0.01, 3.0, 1e-4, Ud);

  std::cout << "Original size " << io_orig.energy_density_clover.size() << " adaptive size " << io_adaptive.energy_density_clover.size() << std::endl;

  //Output values for plotting
  {
    std::ofstream out("t2E_orig.dat");
    out.precision(16);
    for(auto const &e: io_orig.energy_density_clover){
      out << e.first << " " << e.second << std::endl;
    }
  }
  {
    std::ofstream out("t2E_adaptive.dat");
    out.precision(16);
    for(auto const &e: io_adaptive.energy_density_clover){
      out << e.first << " " << e.second << std::endl;
    }
  }

  for(int i=0;i<io_adaptive.energy_density_clover.size();i++){
    RealD t = io_adaptive.energy_density_clover[i].first;
    RealD v_adaptive = io_adaptive.energy_density_clover[i].second;
    RealD v_orig = interpolate(t,  io_orig.energy_density_clover );
    std::cout << t << " orig: " << v_orig << " adaptive: " << v_adaptive << " reldiff: " << (v_adaptive-v_orig)/v_orig << std::endl;
  }

  std::cout << GridLogMessage << " Done" << std::endl;
  Grid_finalize();
  return 0;
}
