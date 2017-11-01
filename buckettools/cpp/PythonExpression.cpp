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
#include "PythonExpression.h"
#include "PythonInstance.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::string &function) : 
                                                dolfin::Expression(), 
                                                pyinst_(function)
{
  assert(pyinst_.number_arguments()==1);
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::size_t &dim, 
                                   const std::string &function) : 
                                            dolfin::Expression(dim), 
                                            pyinst_(function)
{
  assert(pyinst_.number_arguments()==1);
}

//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::vector< std::size_t > 
                                                      &value_shape, 
                                   const std::string &function) : 
                                     dolfin::Expression(value_shape), 
                                     pyinst_(function)
{
  assert(pyinst_.number_arguments()==1);
}

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::string &function, const double_ptr time) : 
                                                dolfin::Expression(), 
                                                pyinst_(function), 
                                                time_(time)
{
  assert((pyinst_.number_arguments()==1)||(pyinst_.number_arguments()==2));
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::size_t &dim, 
                                   const std::string &function,
                                   const double_ptr time) : 
                                            dolfin::Expression(dim), 
                                            pyinst_(function),
                                            time_(time)
{
  assert((pyinst_.number_arguments()==1)||(pyinst_.number_arguments()==2));
}

//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
PythonExpression::PythonExpression(const std::vector< std::size_t > 
                                                      &value_shape, 
                                   const std::string &function,
                                   const double_ptr time) : 
                                     dolfin::Expression(value_shape), 
                                     pyinst_(function),
                                     time_(time)
{
  assert((pyinst_.number_arguments()==1)||(pyinst_.number_arguments()==2));
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
  PyObject *pArgs, *pPos, *pT, *px, *pxx, *pResult;
  std::size_t meshdim;
  
  meshdim = x.size();                                                // the coordinates dimension
  pPos=PyTuple_New(meshdim);                                         // prepare a python tuple for the coordinates
  
  int nargs = pyinst_.number_arguments();
  pArgs = PyTuple_New(nargs);                                        // set up the input arguments tuple
  PyTuple_SetItem(pArgs, 0, pPos);                                   // set the input argument to the coordinates
  if (nargs==2)
  {
    pT=PyFloat_FromDouble(*time_);
    PyTuple_SetItem(pArgs, 1, pT);
  }
  
  if (PyErr_Occurred()){                                             // error check - in setting arguments
    pyinst_.print_error();
    tf_err("In PythonExpression::eval setting pResult.", "Python error occurred.");
  }
  
  for (uint i = 0; i<meshdim; i++)                                   // loop over the coordinate dimensions
  {
    px = PyFloat_FromDouble(x[i]);                                   // convert coordinate to python float
    PyTuple_SetItem(pPos, i, px);                                    // set it in the pPos tuple
  }
  
  pResult = pyinst_.call(pArgs);                                     // call the python function (through the bucket python instance)
  
  if (PyErr_Occurred()){                                             // error check - in running user defined function
    pyinst_.print_error();
    tf_err("In PythonExpression::eval evaluating pResult.", "Python error occurred.");
  }
    
  if (PySequence_Check(pResult))                                     // is the result a sequence
  {                                                                  // yes, ...
    const std::size_t dim0 = PySequence_Length(pResult);
    assert(dim0 <= values.size());
    for (std::size_t i = 0; i < dim0; i++)                           // loop over the value dimension
    {
      px = PySequence_GetItem(pResult, i);                           // get the item from the python sequence
      if (PySequence_Check(px))
      {
        const std::size_t dim1 = PySequence_Length(px);
        assert(dim0*dim1 == values.size());
        for (std::size_t j = 0; j < dim1; j++)
        {
          pxx = PySequence_GetItem(px, j);
          values[i*dim1 + j] = PyFloat_AsDouble(pxx);

          if (PyErr_Occurred()) {
            pyinst_.print_error();
            tf_err("In PythonExpression::eval evaluating tensor values.", "Python error occurred.");
          }

          Py_DECREF(pxx);
        }
      }
      else
      {
        values[i] = PyFloat_AsDouble(px);                            // convert it to a float
        
        if (PyErr_Occurred()){                                       // check for errors in conversion
          pyinst_.print_error();
          tf_err("In PythonExpression::eval evaluating vector values.", "Python error occurred.");
        }
      }
      
      Py_DECREF(px);
    }
  }
  else
  {                                                                  // not a sequence
    values[0] = PyFloat_AsDouble(pResult);                           // just convert a single value
    
    if (PyErr_Occurred()){                                           // check for errors in conversion
      pyinst_.print_error();
      tf_err("In PythonExpression::eval evaluating scalar values.", "Python error occurred.");
    }
  }
  
  Py_DECREF(pResult);                                                // destroy the result python object
  
  Py_DECREF(pArgs);                                                  // destroy the input arugments object
  
}

//*******************************************************************|************************************************************//
// return if this expression is time dependent or not
//*******************************************************************|************************************************************//
const bool PythonExpression::time_dependent() const
{
  return (pyinst_.number_arguments()==2);
}

