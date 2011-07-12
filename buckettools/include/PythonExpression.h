
#ifndef __PYTHON_EXPRESSION_H
#define __PYTHON_EXPRESSION_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"

namespace buckettools
{

  // The PythonExpression class describes a derived dolfin Expression class that overloads
  // the eval function using python
  class PythonExpression : public dolfin::Expression
  {
  // only accessible within this class
  private:
    
    // An instance of python (contains and wraps all the useful python information)
    const PythonInstance pyinst_;

  // accessible to everyone
  public:

    // No default constructor as we require a function

    // Specific constructor
    PythonExpression(const std::string &function);
    
    // Specific constructor (vector)
    PythonExpression(const uint &dim, const std::string &function);
    
    // Specific constructor (tensor)
    PythonExpression(const uint &dim0, const uint &dim1, const std::string &function);
    
    // Specific constructor (alternative tensor)
    PythonExpression(const std::vector<uint> &value_shape, const std::string &function);
    
    // no copy constructor for now
    
    // Default constructor (virtual so calls DOLFIN expression destructor too)
    virtual ~PythonExpression();
    
    // overload eval function to call python instance
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const;
    
  };

}
#endif
