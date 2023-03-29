#pragma once

#include "defines.h"
#include <sys/sysinfo.h>
#include <errno.h>
#include <execinfo.h>

namespace GridMeas{
  using namespace Grid;


  inline double byte_to_MB(const size_t b){
    return double(b)/1024./1024.;
  }

  
  //Print memory information
  void printMem(const std::string &reason = "", bool barrier = true){
    if(barrier) MPI_Barrier(MPI_COMM_WORLD);
    std::cout << GridLogMessage << "====================================" << std::endl;
    std::cout << "printMem called";
    if(reason.size()) std::cout << " with reason: " << reason;
    std::cout << std::endl;
  
    /* unsigned long totalram;  /\* Total usable main memory size *\/ */
    /* unsigned long freeram;   /\* Available memory size *\/ */
    /* unsigned long sharedram; /\* Amount of shared memory *\/ */
    /* unsigned long bufferram; /\* Memory used by buffers *\/ */
    /* unsigned long totalswap; /\* Total swap space size *\/ */
    /* unsigned long freeswap;  /\* swap space still available *\/ */
    /* unsigned short procs;    /\* Number of current processes *\/ */
    /* unsigned long totalhigh; /\* Total high memory size *\/ */
    /* unsigned long freehigh;  /\* Available high memory size *\/ */
    /* unsigned int mem_unit;   /\* Memory unit size in bytes *\/ */

    struct sysinfo myinfo;
    sysinfo(&myinfo);
    double total_mem = myinfo.mem_unit * myinfo.totalram;
    total_mem /= (1024.*1024.);
    double free_mem = myinfo.mem_unit * myinfo.freeram;
    free_mem /= (1024.*1024.);
  
    std::cout << "Memory: total: " << total_mem << " MB, avail: " << free_mem << " MB, used " << total_mem-free_mem << " MB" << std::endl;

#ifdef GRID_CUDA
    size_t gpu_free, gpu_tot;
    cudaError_t err = cudaMemGetInfo(&gpu_free, &gpu_tot);
    if( err != cudaSuccess ) {
      std::cerr << "printMem: cudaMemGetInfo failed: " <<cudaGetErrorString(err)<< std::endl;
      assert(0);
    }
    std::cout << "GPU memory: total: " << byte_to_MB(gpu_tot) << " MB, avail: " << byte_to_MB(gpu_free) << " MB, used " << byte_to_MB(gpu_tot-gpu_free) << " MB" << std::endl; 
#endif

    profilerDebugPrint;

#ifndef GRID_UVM
    MemoryManager::Print();
#endif
    
    std::cout << GridLogMessage << "====================================" << std::endl;
    std::cout.flush();
    if(barrier) MPI_Barrier(MPI_COMM_WORLD);
  }

  //On some devices, stack memory allocations on the device are not freed up as a performance optimization. However, for some kernels that use large amounts of stack memory, this can
  //mean tens of GB that become unusable. This command will reset the cache memory back to the default
  inline void resetDeviceStackMemory(){
#ifdef GRID_CUDA
    cudaDeviceSetLimit(cudaLimitStackSize, 0);
#endif
  }

  
}
