
#include "PythonExpression.h"
#include "PythonInstance.h"
#include <dolfin.h>
#include "Python.h"
#include <string>

using namespace buckettools;

PythonExpression::PythonExpression(std::string function) : dolfin::Expression(), pyinst_(function)
{
  // Do nothing
}
PythonExpression::PythonExpression(uint dim, std::string function) : dolfin::Expression(dim), pyinst_(function)
{
  // Do nothing
}
PythonExpression::PythonExpression(uint dim0, uint dim1, std::string function) : dolfin::Expression(dim0, dim1), pyinst_(function)
{
  // Do nothing
}
PythonExpression::PythonExpression(std::vector<uint> value_shape, std::string function) : dolfin::Expression(value_shape), pyinst_(function)
{
  // Do nothing
}
PythonExpression::~PythonExpression()
{
  // should all be automatic from pyinst_ destructor
}

void PythonExpression::eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
{
  PyObject *pArgs, *pPos, *px, *pResult;
  dolfin::uint meshdim, valdim;
  
  valdim = values.size();
  
  meshdim = x.size();
  pPos=PyTuple_New(meshdim);
  
  pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pPos);
  
  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    dolfin::error("In PythonExpression::eval setting pArgs");
  }
  
  for (dolfin::uint i = 0; i<meshdim; i++)
  {
    px = PyFloat_FromDouble(x[i]);
    PyTuple_SetItem(pPos, i, px);
  }
  
  pResult = pyinst_.call(pArgs);
  
  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    dolfin::error("In PythonExpression::eval evaluating pResult");
  }
    
  if (PySequence_Check(pResult))
  {
    for (dolfin::uint i = 0; i<valdim; i++)
    {
      px = PySequence_GetItem(pResult, i);
      values[i] = PyFloat_AsDouble(px);
      
      // Check for errors in executing user code.
      if (PyErr_Occurred()){
        PyErr_Print();
        dolfin::error("In PythonExpression::eval evaluating values");
      }
      
      Py_DECREF(px);
    }
  }
  else
  {
    values[0] = PyFloat_AsDouble(pResult);
    
    // Check for errors in executing user code.
    if (PyErr_Occurred()){
      PyErr_Print();
      dolfin::error("In PythonExpression::eval evaluating values");
    }
  }
  
  Py_DECREF(pResult);
  
  Py_DECREF(pArgs);
  
}
