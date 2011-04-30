
#ifndef __DETECTORS_H
#define __DETECTORS_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"

namespace buckettools
{
  
  typedef boost::shared_ptr< dolfin::Array<double> > Array_double_ptr;
  
  class Detectors
  {
  protected:
    
    std::string  name_;
    dolfin::uint number_detectors_, 
                 meshdim_;
    std::vector< Array_double_ptr > positions_;
    
    void clean_();
    
  public:
    // no default constructor for now
    
    // Constructor - takes number of detectors and a python function string
    Detectors();
    
    Detectors(dolfin::uint number_detectors, 
              dolfin::uint meshdim, 
              std::string  name);
    
    // no copy constructor for now
    
    // Destructor
    ~Detectors();
    
    void eval(std::vector< Array_double_ptr > &values,
              dolfin::GenericFunction         &function);
    
    std::string str() const;
    
    std::vector< Array_double_ptr >::iterator begin();

    std::vector< Array_double_ptr >::const_iterator begin() const;

    std::vector< Array_double_ptr >::iterator end();

    std::vector< Array_double_ptr >::const_iterator end() const;
    
    std::string name() const;
    
    uint dim() const;
    
    uint size() const;

  };
  
  typedef boost::shared_ptr< Detectors > Detectors_ptr;
  
}

#endif