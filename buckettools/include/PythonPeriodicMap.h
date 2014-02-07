
#ifndef __PYTHON_PERIODIC_MAP_H
#define __PYTHON_PERIODIC_MAP_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // PythonPeriodicMap class:
  //
  // The PythonPeriodicMap class describes a derived dolfin SubDomain class that overloads
  // the map function using python
  //*****************************************************************|************************************************************//
  class PythonPeriodicMap : public dolfin::SubDomain
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    PythonPeriodicMap(const std::string &function);                  // specific constructor
    
    virtual ~PythonPeriodicMap();                                    // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void map(const dolfin::Array<double>& x, dolfin::Array<double>& y) const;        // map slave position x to master position y
    

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//
    
    const PythonInstance pyinst_;                                    // a python instance (wrapping useful python information)

  };

  typedef boost::shared_ptr< PythonPeriodicMap > PythonPeriodicMap_ptr;// define a (boost shared) pointer type for the python
                                                                     // periodic map

}
#endif
