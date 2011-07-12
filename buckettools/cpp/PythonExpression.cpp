
#include "PythonExpression.h"
#include "PythonInstance.h"
#include <dolfin.h>
#include "Python.h"
#include <string>

using namespace buckettools;

// Specific constructor
PythonExpression::PythonExpression(const std::string &function) : dolfin::Expression(), pyinst_(function)
{
  // Do nothing
}

// Specific constructor (vector)
PythonExpression::PythonExpression(const uint &dim, const std::string &function) : dolfin::Expression(dim), pyinst_(function)
{
  // Do nothing
}

// Specific constructor (tensor)
PythonExpression::PythonExpression(const uint &dim0, const uint &dim1, const std::string &function) : dolfin::Expression(dim0, dim1), pyinst_(function)
{
  // Do nothing
}

// Specific constructor (alternative tensor)
PythonExpression::PythonExpression(const std::vector<uint> &value_shape, const std::string &function) : dolfin::Expression(value_shape), pyinst_(function)
{
  // Do nothing
}

// Default destructor
PythonExpression::~PythonExpression()
{
  // should all be automatic from pyinst_ destructor
}

// Overload dolfin eval
void PythonExpression::eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
{
  PyObject *pArgs, *pPos, *px, *pResult;
  dolfin::uint meshdim, valdim;
  
  // the size of the value space (sadly doesn't tell us about shape so tensors not quite supported yet)
  valdim = values.size();
  
  // the mesh dimension
  meshdim = x.size();
  pPos=PyTuple_New(meshdim);
  
  // set up the arguments to the val function (only x so far)
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
    
  // Is the result a sequence (assumed vector for now)?
  if (PySequence_Check(pResult))
  {
    // yes, then loop over valdim filling in values
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
    // no, just get the scalar value back then
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

