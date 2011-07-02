
#ifndef __GENERICDETECTORS_H
#define __GENERICDETECTORS_H

#include <dolfin.h>

namespace buckettools
{
  
  typedef boost::shared_ptr< dolfin::Array<double> > Array_double_ptr;
  
  // GenericDetectors is a generic class that contains base information for point and python detectors
  class GenericDetectors
  {
  protected:
    
    // The name of this detector set
    std::string                     name_;
    // The number of detectors and the mesh dimension
    dolfin::uint                    number_detectors_, 
                                    meshdim_;
    // A vector of arrays describing the locations of the positions
    std::vector< Array_double_ptr > positions_;
    
    // Clean the detector class
    void clean_();
    
  public:
    // Default constructor
    GenericDetectors();
    
    // Specific constructor - takes number of detectors, the mesh dimension and the detector's name
    GenericDetectors(dolfin::uint number_detectors, 
              dolfin::uint meshdim, 
              std::string  name);
    
    // No copy constructor yet
    
    // Destructor
    ~GenericDetectors();
    
    // The base implementation of eval, which returns the values of a function at the locations described by positions_
    void eval(std::vector< Array_double_ptr > &values,
              dolfin::GenericFunction         &function);
    
    // Return the name_
    std::string name() const;
    
    // Return the mesh dimension
    uint dim() const;
    
    // Return how many detectors are in this set
    uint size() const;

    // Pretty print output of the detector positions_
    std::string str() const;
    
    // Return begin iterator for positions_
    std::vector< Array_double_ptr >::iterator begin();

    // Return const begin iterator for positions_
    std::vector< Array_double_ptr >::const_iterator begin() const;

    // Return end iterator for positions_
    std::vector< Array_double_ptr >::iterator end();

    // Return const end iterator for positions_
    std::vector< Array_double_ptr >::const_iterator end() const;
    
  };
  
  // Define a boost shared_ptr for the base GenericDetectors class
  typedef boost::shared_ptr< GenericDetectors > GenericDetectors_ptr;
  
}

#endif
