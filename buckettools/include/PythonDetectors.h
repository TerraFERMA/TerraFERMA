
#ifndef __PYTHON_DETECTORS_H
#define __PYTHON_DETECTORS_H

#include "Detectors.h"
#include <dolfin.h>
#include "PythonInstance.h"

namespace buckettools
{
  
  class PythonDetectors : public Detectors
  {
  private:
    
    PythonInstance pyinst_;
    
    void init_();
    
  public:
    // no default constructor for now
    
    // Constructor - takes number of detectors and a python function string
    PythonDetectors(dolfin::uint number_detectors, dolfin::uint meshdim, std::string function, std::string name);
    
    // no copy constructor for now
    
    // Destructor
    virtual ~PythonDetectors();
    
    virtual std::string str() const;
    
  };
  
  typedef boost::shared_ptr< PythonDetectors > PythonDetectors_ptr;
  
}

#endif
