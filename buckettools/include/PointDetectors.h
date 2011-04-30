
#ifndef __POINT_DETECTORS_H
#define __POINT_DETECTORS_H

#include <dolfin.h>
#include "Detectors.h"

namespace buckettools
{
  
  class PointDetectors : public Detectors
  {
  private:
    
    void init_(Array_double_ptr point);
    
    void init_(std::vector<double> *point);
    
  public:
    PointDetectors();
    
    // Constructor
    PointDetectors(Array_double_ptr point, std::string name);
    
    // Constructor
    PointDetectors(std::vector<double> *point, std::string name);
    
    // no copy constructor for now
    
    // Destructor
    virtual ~PointDetectors();
    
  };
  
  typedef boost::shared_ptr< PointDetectors > PointDetectors_ptr;
  
}

#endif