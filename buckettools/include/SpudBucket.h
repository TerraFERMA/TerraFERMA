
#ifndef __SPUD_BUCKET_H
#define __SPUD_BUCKET_H

#include "Bucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  class SpudBucket : public Bucket
  {
  private:
    
    void meshes_fill_(const std::string &optionpath);
    
    void system_fill_(const std::string &optionpath);
    
    void detectors_fill_();
    
    void bc_fill_(const std::string bcpath, 
                  const int funci,
                  const int dimension,
                  FunctionSpace_ptr subsysspace,
                  const std::string sysname,
                  const std::string funcname);
    
    GenericFunction_ptr init_exp_(const std::string path, const int dimension);
    
  public:
    
    SpudBucket()
    { SpudBucket("uninitialised_name", "uninitialised_path"); }
    
    SpudBucket(std::string name, std::string option_path);
    
    virtual ~SpudBucket();
    
    void fill();
    
  };
}
#endif
