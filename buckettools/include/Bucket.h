
#ifndef __BUCKET_H
#define __BUCKET_H

#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>

namespace buckettools
{
  
  typedef boost::shared_ptr< dolfin::GenericFunction > GenericFunction_ptr;
  
  class Bucket
  {
  private:
    
    std::vector< Detectors_ptr > detectors_;
    
    std::vector< GenericFunction_ptr > functions_;
    
    void clean_();
    
  public:
    
    Bucket();
    
    ~Bucket();
    
    void register_detector(Detectors_ptr detector);
    
    void register_function(GenericFunction_ptr function);
    
    std::vector< Detectors_ptr >::iterator detectors_begin();
    
    std::vector< Detectors_ptr >::const_iterator detectors_begin() const;
    
    std::vector< Detectors_ptr >::iterator detectors_end();
    
    std::vector< Detectors_ptr >::const_iterator detectors_end() const;
    
    std::vector< GenericFunction_ptr >::iterator functions_begin();
    
    std::vector< GenericFunction_ptr >::const_iterator functions_begin() const;
    
    std::vector< GenericFunction_ptr >::iterator functions_end();
    
    std::vector< GenericFunction_ptr >::const_iterator functions_end() const;
    
  };
}
#endif