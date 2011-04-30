
#ifndef __PYTHON_EXPRESSION_H
#define __PYTHON_EXPRESSION_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"

namespace buckettools
{

  class PythonExpression : public dolfin::Expression
  {
  private:
    
    const PythonInstance pyinst_;
    
  public:
    PythonExpression(std::string function);
    
    PythonExpression(uint dim, std::string function);
    
    PythonExpression(uint dim0, uint dim1, std::string function);
    
    PythonExpression(std::vector<uint> value_shape, std::string function);
    
    // no copy constructor for now
    
    virtual ~PythonExpression();
    
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const;
    
  };

}
#endif