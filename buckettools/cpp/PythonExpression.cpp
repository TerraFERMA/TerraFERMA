
#include "PythonExpression.h"
#include "PythonInstance.h"
#include <dolfin.h>
#include "Python.h"
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::string &function) : 
                                                dolfin::Expression(), 
                                                pyinst_(function)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const uint &dim, 
                                   const std::string &function) : 
                                            dolfin::Expression(dim), 
                                            pyinst_(function)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const uint &dim0, 
                                   const uint &dim1, 
                                   const std::string &function) : 
                                      dolfin::Expression(dim0, dim1),
                                      pyinst_(function)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (alternate tensor)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::vector<uint> 
                                                      &value_shape, 
                                   const std::string &function) : 
                                     dolfin::Expression(value_shape), 
                                     pyinst_(function)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PythonExpression::~PythonExpression()
{
                                                                     // should all be automatic from pyinst_ destructor
}

//*******************************************************************|************************************************************//
// overload dolfin eval
//*******************************************************************|************************************************************//
void PythonExpression::eval(dolfin::Array<double>& values, 
                            const dolfin::Array<double>& x) const
{
  PyObject *pArgs, *pPos, *px, *pResult;
  dolfin::uint meshdim, valdim;
  
  valdim = values.size();                                            // the size of the value space 
                                                                     // (FIXME: doesn't tell us about shape so tensors not supported)
  
  meshdim = x.size();                                                // the coordinates dimension
  pPos=PyTuple_New(meshdim);                                         // prepare a python tuple for the coordinates
  
  pArgs = PyTuple_New(1);                                            // set up the input arguments tuple (just 1 argument now - can
                                                                     // we tell this from what the user writes and decide how many 
                                                                     // this should be?)
  PyTuple_SetItem(pArgs, 0, pPos);                                   // set the input argument to the coordinates
  
  if (PyErr_Occurred()){                                             // error check - in setting arguments
    PyErr_Print();
    dolfin::error("In PythonExpression::eval setting pArgs");
  }
  
  for (uint i = 0; i<meshdim; i++)                                   // loop over the coordinate dimensions
  {
    px = PyFloat_FromDouble(x[i]);                                   // convert coordinate to python float
    PyTuple_SetItem(pPos, i, px);                                    // set it in the pPos tuple
  }
  
  pResult = pyinst_.call(pArgs);                                     // call the python function (through the bucket python instance)
  
  if (PyErr_Occurred()){                                             // error check - in running user defined function
    PyErr_Print();
    dolfin::error("In PythonExpression::eval evaluating pResult");
  }
    
  if (PySequence_Check(pResult))                                     // is the result a sequence (FIXME: assumed vector, not tensor)
  {                                                                  // yes, ...
    for (dolfin::uint i = 0; i<valdim; i++)                          // loop over the value dimension
    {
      px = PySequence_GetItem(pResult, i);                           // get the item from the python sequence
      values[i] = PyFloat_AsDouble(px);                              // convert it to a float
      
      if (PyErr_Occurred()){                                         // check for errors in conversion
        PyErr_Print();
        dolfin::error("In PythonExpression::eval evaluating values");
      }
      
      Py_DECREF(px);
    }
  }
  else
  {                                                                  // not a sequence
    values[0] = PyFloat_AsDouble(pResult);                           // just convert a single value
    
    if (PyErr_Occurred()){                                           // check for errors in conversion
      PyErr_Print();
      dolfin::error("In PythonExpression::eval evaluating values");
    }
  }
  
  Py_DECREF(pResult);                                                // destroy the result python object
  
  Py_DECREF(pArgs);                                                  // destroy the input arugments object
  
}

