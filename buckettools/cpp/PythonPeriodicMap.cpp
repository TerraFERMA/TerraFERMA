
#include "PythonPeriodicMap.h"
#include "PythonInstance.h"
#include <dolfin.h>
#include "Python.h"
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PythonPeriodicMap::PythonPeriodicMap(const std::string &function) : 
                                                dolfin::SubDomain(), 
                                                pyinst_(function)
{
  assert(pyinst_.number_arguments()==1);
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PythonPeriodicMap::~PythonPeriodicMap()
{
                                                                     // should all be automatic from pyinst_ destructor
}

//*******************************************************************|************************************************************//
// overload dolfin map
//*******************************************************************|************************************************************//
void PythonPeriodicMap::map(const dolfin::Array<double>& x, 
                            dolfin::Array<double>& y) const
{
  PyObject *pArgs, *pPos, *pT, *px, *pResult;
  std::size_t meshdim, valdim;
  
  valdim = y.size();
  meshdim = x.size();                                                // the coordinates dimension
  assert(valdim==meshdim);

  pPos=PyTuple_New(meshdim);                                         // prepare a python tuple for the coordinates
  
  pArgs = PyTuple_New(1);                                            // set up the input arguments tuple
  PyTuple_SetItem(pArgs, 0, pPos);                                   // set the input argument to the coordinates
  
  if (PyErr_Occurred()){                                             // error check - in setting arguments
    PyErr_Print();
    dolfin::error("In PythonPeriodicMap::map setting pArgs");
  }
  
  for (uint i = 0; i<meshdim; i++)                                   // loop over the coordinate dimensions
  {
    px = PyFloat_FromDouble(x[i]);                                   // convert coordinate to python float
    PyTuple_SetItem(pPos, i, px);                                    // set it in the pPos tuple
  }
  
  pResult = pyinst_.call(pArgs);                                     // call the python function (through the bucket python instance)
  
  if (PyErr_Occurred()){                                             // error check - in running user defined function
    PyErr_Print();
    dolfin::error("In PythonPeriodicMap::map evaluating pResult");
  }
    
  if (PySequence_Check(pResult))                                     // is the result a sequence
  {                                                                  // yes, ...
    for (std::size_t i = 0; i<valdim; i++)                          // loop over the value dimension
    {
      px = PySequence_GetItem(pResult, i);                           // get the item from the python sequence
      y[i] = PyFloat_AsDouble(px);                                   // convert it to a float
      
      if (PyErr_Occurred()){                                         // check for errors in conversion
        PyErr_Print();
        dolfin::error("In PythonPeriodicMap::map evaluating values");
      }
      
      Py_DECREF(px);
    }
  }
  else
  {                                                                  // not a sequence
    y[0] = PyFloat_AsDouble(pResult);                                // just convert a single value
    
    if (PyErr_Occurred()){                                           // check for errors in conversion
      PyErr_Print();
      dolfin::error("In PythonPeriodicMap::map evaluating values");
    }
  }
  
  Py_DECREF(pResult);                                                // destroy the result python object
  
  Py_DECREF(pArgs);                                                  // destroy the input arugments object
  
}

