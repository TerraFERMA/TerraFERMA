
#ifndef __DIAGNOSTICS_FILE_H
#define __DIAGNOSTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "StatFile.h"
#include "Bucket.h"

namespace buckettools
{

  class DiagnosticsFile : public StatFile
  {
  public:
    
    DiagnosticsFile(const std::string name);
 
    ~DiagnosticsFile();
    
    void write_header(const Bucket &bucket, 
                      const bool &timestepping);

    void write_data(Bucket &bucket);
    
    void write_data(const uint   &timestep,
                    const double &elapsedtime, 
                    const double &dt, 
                    Bucket       &bucket);
    
  private:

    void header_bucket_(const Bucket &bucket,
                        uint &column);

    void header_field_(FunctionBucket_const_it f_begin, 
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_coeff_(FunctionBucket_const_it f_begin, 
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_functional_(FunctionBucket_const_it f_it,
                            Form_const_it s_begin,
                            Form_const_it s_end, 
                            uint &column);

    void data_timestep_(const uint   &timestep,
                        const double &elapsedtime, 
                        const double &dt);
    
    void data_bucket_(Bucket &bucket);

    void data_field_(FunctionBucket_const_it f_begin, 
                     FunctionBucket_const_it f_end);

    void data_coeff_(FunctionBucket_const_it f_begin, 
                     FunctionBucket_const_it f_end);

    void data_functional_(Form_const_it s_begin,
                          Form_const_it s_end);

  };
  
}
#endif
