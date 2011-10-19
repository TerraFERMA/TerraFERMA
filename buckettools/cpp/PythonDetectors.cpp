
#include "PythonDetectors.h"
#include "GenericDetectors.h"
#include <dolfin.h>
#include "Python.h"
#include <string>
#include "PythonInstance.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PythonDetectors::PythonDetectors(const uint &number_detectors, 
                                 const uint &meshdim, 
                                 const std::string &function, 
                                 const std::string &name) : 
                  GenericDetectors(number_detectors, meshdim, name), 
                  pyinst_(function)
{
  init_();                                                           // initialize
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PythonDetectors::~PythonDetectors()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// intialize a python detector from the pyinst_
//*******************************************************************|************************************************************//
void PythonDetectors::init_()
{

  assert(pyinst_.number_arguments()==0);

  if(!positions_.empty())
  {
    dolfin::error("In PythonDetectors::init_ intializing already initialized detectors.");
  }
                                                                     // python objects:
  PyObject         *pArgs,                                           // input arguments
                   *px,                                              // coordinate
                   *pResult,                                         // result
                   *pResultItem;                                     // result component
  Array_double_ptr point;                                            // set up a detector location to be reset
  
  pArgs = PyTuple_New(0);                                            // set up the input arguments (contain nothing!)
  
  pResult = pyinst_.call(pArgs);                                     // run the python function through the python instance
  
  if (PyErr_Occurred()){                                             // check for errors evaluating user code
    PyErr_Print();
    dolfin::error("In PythonDetectors::init_ evaluating pResult");
  }
    
  for (dolfin::uint i = 0; i<number_detectors_; i++)                 // loop over the array of detectors
  {
    pResultItem = PySequence_GetItem(pResult, i);                    // get an item out of the sequence of the results
    
    point.reset(new dolfin::Array<double>(meshdim_));                // set the item as a new dolfin array
    
    for (dolfin::uint j = 0; j<meshdim_; j++)                        // loop over the coordinate dimension
    {
      px = PySequence_GetItem(pResultItem, j);                       // get an item from the coordinate array
      (*point)[j] = PyFloat_AsDouble(px);                            // convert from python to double
      
      if (PyErr_Occurred()){                                         // check for errors in conversion
        PyErr_Print();
        dolfin::error("In PythonDetectors::init_ evaluating values");
      }
      
      Py_DECREF(px);                                                 // deallocate python object
    }
    
    positions_.push_back(point);                                     // save the point
    
    Py_DECREF(pResultItem);                                          // deallocate python object
    
  }
  
  Py_DECREF(pResult);                                                // deallocate python object
  Py_DECREF(pArgs);                                                  // deallocate python object
}

//*******************************************************************|************************************************************//
// return a string describing the positions of the detectors and the python function
//*******************************************************************|************************************************************//
const std::string PythonDetectors::str() const
{
  std::stringstream s;
  
  s << pyinst_.function() << std::endl;
  
  s << GenericDetectors::str();
  
  return s.str();
}

