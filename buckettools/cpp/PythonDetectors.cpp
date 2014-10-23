// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#include "Python.h"
#include "PythonDetectors.h"
#include "GenericDetectors.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include "PythonInstance.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PythonDetectors::PythonDetectors(const uint &meshdim, 
                                 const std::string &function, 
                                 const std::string &name) : 
                  GenericDetectors(-1, meshdim, name),               // don't know size yet
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
    tf_err("Intializing already initialized detectors.", "positions_.empty() not empty.");
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
    tf_err("In PythonDetectors::init_ evaluating pResult.", "Python error occurred.");
  }

  number_detectors_ = PyObject_Length(pResult);                      // find out how many detectors we have
  assert(size()!=-1);
  assert(size()>0);
    
  if (PyErr_Occurred()){                                             // check for errors evaluating user code
    PyErr_Print();
    tf_err("In PythonDetectors::init_ evaluating length of pResult.", "Python error occurred.");
  }
    
  for (std::size_t i = 0; i<size(); i++)                             // loop over the array of detectors
  {
    pResultItem = PySequence_GetItem(pResult, i);                    // get an item out of the sequence of the results
    
    point.reset(new dolfin::Array<double>(meshdim_));                // set the item as a new dolfin array

    assert(PyObject_Length(pResultItem)==meshdim_);                  // check this item in the list is meshdim_ long
    
    for (std::size_t j = 0; j<meshdim_; j++)                         // loop over the coordinate dimension
    {
      px = PySequence_GetItem(pResultItem, j);                       // get an item from the coordinate array
      (*point)[j] = PyFloat_AsDouble(px);                            // convert from python to double
      
      if (PyErr_Occurred()){                                         // check for errors in conversion
        PyErr_Print();
        tf_err("In PythonDetectors::init_ evaluating values.", "Python error occurred.");
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

