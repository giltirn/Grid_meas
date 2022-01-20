#pragma once

#include "field_array_io.h"

namespace GridMeas{
  using namespace Grid;

  //Binary write of real array
  //Use only one a single node
  inline void binaryWriteArray(const std::string &file, const std::vector<RealD> &data){
    std::string format("IEEE64BIG");    
    std::vector<RealD> cp(data);
    assert(sizeof(RealD) == 8);
    size_t size = data.size();
    BinaryIO::htobe64_v(cp.data(), size*8);
    
    uint32_t checksum = 0;
    uint32_t const* idata = (uint32_t const*)cp.data();
    
    for(size_t i=0;i<2*size;i++)
      mergeChecksum(checksum, idata[i]);
    
    NerscIO::truncate(file);

    size_t offset = writeFieldArrayHeader(size, checksum, format, file);
    
    std::ofstream fout(file,std::ios::out|std::ios::in);
    fout.seekp(offset,std::ios::beg);

    fout.write((char const*)cp.data(), size* 8);
  }

  //Binary read of real array
  inline void binaryReadArray(std::vector<RealD> &data, const std::string &file){
    uint32_t size, hdr_checksum;
    std::string format;
    size_t offset = readFieldArrayHeader(size, hdr_checksum, format, file);

    assert(format == "IEEE64BIG");
    data.resize(size);

    std::ifstream in(file);
    in.seekg(offset);
    
    in.read((char*)data.data(), size * 8);
    
    uint32_t checksum = 0;
    uint32_t const* idata = (uint32_t const*)data.data();
    
    for(size_t i=0;i<2*size;i++)
      mergeChecksum(checksum, idata[i]);

    assert(checksum == hdr_checksum);

    BinaryIO::be64toh_v((void*)data.data(), size * 8);
  }

  //Write a real number to plain text
  void asciiWriteValue(const RealD data, const std::string &stub, const int traj){
    std::string filename = stub + "." + std::to_string(traj);
    std::ofstream of(filename);
    of.precision(17);
    of << data << std::endl;
  }

  //Write a real array to plain text on the same line
  void asciiWriteArray(const std::vector<RealD> &data, const std::string &stub, const int traj){
    std::string filename = stub + "." + std::to_string(traj);
    std::ofstream of(filename);
    of.precision(17);

    for(int t=0;t<data.size();t++){
      of << data[t] << (t < data.size() - 1 ? " " : "");
    }
    of << std::endl;
  }

  //Write a real array of x,y pairs to plain text on different lines
  void asciiWriteArray(const std::vector<std::pair<RealD,RealD> > &data, const std::string &stub, const int traj){
    std::string filename = stub + "." + std::to_string(traj);
    std::ofstream of(filename);
    of.precision(17);

    for(int t=0;t<data.size();t++){
      of << data[t].first << " " << data[t].second << std::endl;      
    }
  }

  //Write a real array-of-arrays to plain text with all data for a given outer index on the same line
  void asciiWriteArray(const std::vector<std::vector<RealD> > &data, const std::string &stub, const int traj){
    std::string filename = stub + "." + std::to_string(traj);
    std::ofstream of(filename);
    of.precision(17);

    for(int i=0;i<data.size();i++){
      for(int t=0;t<data[i].size();t++){
	of << data[i][t] << (t < data.size() - 1 ? " " : "");
      }
      of << std::endl;
    }
  }

  bool fileExists(const std::string &fn){
    std::ifstream f(fn);
    return f.good();
  }
    
};
