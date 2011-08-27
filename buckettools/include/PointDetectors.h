
#ifndef __POINT_DETECTORS_H
#define __POINT_DETECTORS_H

#include <dolfin.h>
#include "GenericDetectors.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // PointDetectors class:
  //
  // PointDetectors is a derived class of GenericDetectors that implements a detector at a single point
  //*****************************************************************|************************************************************//
  class PointDetectors : public GenericDetectors
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    PointDetectors();                                                // generic constructor
    
    PointDetectors(const Array_double_ptr point, 
                                    const std::string &name);        // specific constructor (dolfin array)
    
    PointDetectors(const std::vector<double> &point, 
                                        const std::string &name);    // specific constructor (std vector)
    
    ~PointDetectors();                                               // default destructor
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_(const Array_double_ptr point);                        // initialize the point detector at a point defined by an array
    
    void init_(const std::vector<double> &point);                    // initialize the point detector at a point defined by a vector
    
  };
  
  typedef boost::shared_ptr< PointDetectors > PointDetectors_ptr;    // define a (boost shared) pointer for this class type
  
}

#endif
