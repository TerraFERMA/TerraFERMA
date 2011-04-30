
#ifndef __DETECTORS_FILE_H
#define __DETECTORS_FILE_H

#include "Bucket.h"
#include <cstdio>
#include <fstream>
#include <string>
#include "StatFile.h"

namespace buckettools
{

  class DetectorsFile : public StatFile
  {
  public:
    
    DetectorsFile(const std::string name);
    
    void write_header(Bucket &bucket,
                      const bool timestepping);
    
    void write_data(Bucket &bucket);
    
    void write_data(const uint   timestep,
                    const double elapsedtime, 
                    const double dt, 
                    Bucket       &bucket);
    
  private:
    
    void header_timestep_(uint &column);
    
    void header_bucket_(Bucket &bucket,
                        uint &column);
    
    void data_timestep_(const uint   timestep,
                        const double elapsedtime, 
                        const double dt);
    
    void data_bucket_(Bucket &bucket);

  };
  
}
#endif