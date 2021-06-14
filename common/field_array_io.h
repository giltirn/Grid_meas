#pragma once

#include<Grid/Grid.h>

namespace GridMeas{
  using namespace Grid;

  inline void mergeChecksum(uint32_t &into, const uint32_t what){
    into ^= what + 0x9e3779b9 + (into<<6) + (into>>2);
  }

  inline uint64_t writeFieldArrayHeader(const uint32_t size, const uint32_t checksum, const std::string &format, const std::string &file){
    std::ofstream fout(file,std::ios::out|std::ios::in);
    fout.seekp(0,std::ios::beg);
    fout << std::setw(10) << size << std::endl;
    fout << std::hex << std::setw(10) << checksum << std::endl;
    fout << format << std::endl;
    return fout.tellp();
  }
 
  inline uint64_t readFieldArrayHeader(uint32_t &size, uint32_t &checksum, std::string &format, const std::string &file){
    std::ifstream fin(file);
    std::string line;
    getline(fin,line);
    {
      std::stringstream ss; ss <<line ; ss >> size;
    }
    getline(fin,line);
    {
      std::stringstream ss; ss <<line ; ss >> std::hex >> checksum;
    }
    getline(fin,format);
    removeWhitespace(format);
      
    return fin.tellg();
  }
 
  template<typename FieldType>
  void writeFieldArray(const std::string &file, const std::vector<FieldType> &data){
    typedef typename FieldType::vector_object vobj;
    typedef typename FieldType::scalar_object sobj;
    GridBase* grid = data[0].Grid(); //assume all fields have the same Grid
    BinarySimpleMunger<sobj, sobj> munge; //straight copy

    //We need a 2-pass header write, first to establish the size, the second pass writes the checksum
    std::string format = getFormatString<typename FieldType::vector_object>();

    uint64_t offset; //leave 64 bits for header
    if ( grid->IsBoss() ) { 
      NerscIO::truncate(file);
      offset = writeFieldArrayHeader(data.size(), 0, format, file);
    }
    grid->Broadcast(0,(void *)&offset,sizeof(offset)); //use as a barrier

    std::cout << "Data offset write " << offset << std::endl;
    std::cout << "Data size write " << data.size() << std::endl;
    uint64_t field_size = uint64_t(grid->gSites()) * sizeof(sobj);
    std::cout << "Field size = " << field_size << " B" << std::endl;

    uint32_t checksum = 0;
    for(int i=0;i<data.size();i++){
      std::cout << "Data field write " << i << " offset " << offset << std::endl;
      uint32_t nersc_csum,scidac_csuma,scidac_csumb;
      BinaryIO::writeLatticeObject<vobj,sobj>(const_cast<FieldType &>(data[i]),file,munge,offset,format,
					      nersc_csum,scidac_csuma,scidac_csumb);
      offset += field_size;
      mergeChecksum(checksum, nersc_csum);
    }
    std::cout << "Write checksum " << checksum << std::endl;

    if ( grid->IsBoss() ) { 
      writeFieldArrayHeader(data.size(), checksum, format, file);
    }
  }


  template<typename FieldType>
  void readFieldArray(std::vector<FieldType> &data, const std::string &file){
    typedef typename FieldType::vector_object vobj;
    typedef typename FieldType::scalar_object sobj;
    assert(data.size() > 0);
    GridBase* grid = data[0].Grid(); //assume all fields have the same Grid
    BinarySimpleUnmunger<sobj, sobj> munge; //straight copy
  
    uint32_t hdr_checksum, hdr_size;
    std::string format;
    uint64_t offset = readFieldArrayHeader(hdr_size, hdr_checksum, format, file);
  
    std::cout << "Data offset read " << offset << std::endl;  
    std::cout << "Data size read " << hdr_size << std::endl;
    assert(data.size() == hdr_size);

    uint64_t field_size = uint64_t(grid->gSites()) * sizeof(sobj);

    uint32_t checksum = 0;

    for(int i=0;i<data.size();i++){
      std::cout << "Data field read " << i << " offset " << offset << std::endl;
      uint32_t nersc_csum,scidac_csuma,scidac_csumb;
      BinaryIO::readLatticeObject<vobj,sobj>(data[i],file,munge,offset,format,
					     nersc_csum,scidac_csuma,scidac_csumb);
      offset += field_size;
      mergeChecksum(checksum, nersc_csum);
    }

    std::cout << "Header checksum " << hdr_checksum << std::endl;    
    std::cout << "Read checksum " << checksum << std::endl;
    

    assert( hdr_checksum == checksum );
  }

};
