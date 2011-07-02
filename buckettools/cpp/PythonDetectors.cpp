
#include "PythonDetectors.h"
#include "GenericDetectors.h"
#include <dolfin.h>
#include "Python.h"
#include <string>
#include "PythonInstance.h"

using namespace buckettools;

PythonDetectors::PythonDetectors(dolfin::uint number_detectors, 
                                 dolfin::uint meshdim, 
                                 std::string function, 
                                 std::string name) : GenericDetectors(number_detectors, meshdim, name), pyinst_(function)
{
  init_();
}

PythonDetectors::~PythonDetectors()
{
  // Do nothing - pyinst_ destructor takes care of it automatically (hopefully) and the base class destructor is virtual
}

void PythonDetectors::init_()
{
  if(!positions_.empty())
  {
    dolfin::error("In PythonDetectors::init_ intializing already initialized detectors.");
  }
  
  PyObject         *pArgs, 
                   *px, 
                   *pResult, 
                   *pResultItem;
  Array_double_ptr point(new dolfin::Array<double>(meshdim_));
  
  pArgs = PyTuple_New(0);
  
  pResult = pyinst_.call(pArgs);
  
  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    dolfin::error("In PythonDetectors::init_ evaluating pResult");
  }
    
  for (dolfin::uint i = 0; i<number_detectors_; i++)
  {
    pResultItem = PySequence_GetItem(pResult, i);
    
    point.reset(new dolfin::Array<double>(meshdim_));
    
    for (dolfin::uint j = 0; j<meshdim_; j++)
    {
      px = PySequence_GetItem(pResultItem, j);
      (*point)[j] = PyFloat_AsDouble(px);
      
      // Check for errors in executing user code.
      if (PyErr_Occurred()){
        PyErr_Print();
        dolfin::error("In PythonDetectors::init_ evaluating values");
      }
      
      Py_DECREF(px);      
    }
    
    positions_.push_back(point);
    
    Py_DECREF(pResultItem);
    
  }
  
  Py_DECREF(pResult);
  
  Py_DECREF(pArgs);  
}

std::string PythonDetectors::str() const
{
  std::stringstream s;
  
  s << pyinst_.function() << std::endl;
  
  s << GenericDetectors::str();
  
  return s.str();
}
