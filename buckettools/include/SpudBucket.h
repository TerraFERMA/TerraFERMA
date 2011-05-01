
#ifndef __SPUD_BUCKET_H
#define __SPUD_BUCKET_H

#include "Bucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  class SpudBucket : public Bucket
  {
  private:
    
    void meshes_fill_();
    
    void detectors_fill_();
    
    void system_fill_(const uint &sysindex);
    
    void bc_fill_(const std::string bcpath, 
                  const int funci,
                  const int dimension,
                  FunctionSpace_ptr subsysspace,
                  const std::string sysname,
                  const std::string funcname);
    
    GenericFunction_ptr init_exp_(const std::string path, const int dimension);
    
  public:
    
    SpudBucket();
    
    virtual ~SpudBucket();
    
    void fill();
    
  };
}
#endif