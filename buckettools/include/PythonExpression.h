
#ifndef __PYTHON_EXPRESSION_H
#define __PYTHON_EXPRESSION_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // PythonExpression class:
  //
  // The PythonExpression class describes a derived dolfin Expression class that overloads
  // the eval function using python
  //*****************************************************************|************************************************************//
  class PythonExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    PythonExpression(const std::string &function);                   // specific constructor (scalar)
    
    PythonExpression(const uint &dim, const std::string &function);  // specific constructor (vector)
    
    PythonExpression(const uint &dim0, const uint &dim1,             // specific constructor (tensor)
                     const std::string &function);
    
    PythonExpression(const std::vector<uint> &value_shape,           // specific constructor (alternate tensor)
                     const std::string &function);
    
    virtual ~PythonExpression();                                     // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate the expression at a given point
              const dolfin::Array<double>& x) const;                 // (but no cell information?)
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//
    
    const PythonInstance pyinst_;                                    // a python instance (wrapping useful python information)

  };

}
#endif
