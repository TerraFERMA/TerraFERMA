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


#include "GlobalPythonInstance.h"
#include "Python.h"
#include "PythonInstance.h"
#include "Logger.h"
#include <string>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PythonInstance::PythonInstance(const std::string &function) : 
                                                  function_(function)
{
  init_();                                                           // initialize
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PythonInstance::~PythonInstance()
{
  clean_();                                                          // finalize
}

//*******************************************************************|************************************************************//
// given a python arguments object call the val function and return the result
//*******************************************************************|************************************************************//
PyObject* PythonInstance::call(PyObject *pArgs) const
{
  return PyObject_CallObject(pFunc_, pArgs);
}

//*******************************************************************|************************************************************//
// print an error message
//*******************************************************************|************************************************************//
void PythonInstance::print_error() const
{
  log(ERROR, "Python computation raised an exception.");
  unsigned int lineno = 0;
  std::istringstream functionss(function_);
  std::stringstream functionstream; functionstream.str("");
  functionstream << std::string(80, '-') << std::endl;
  std::string line;
  while (std::getline(functionss, line)) 
  {
    functionstream << std::setw(4) << ++lineno << "  " << line << std::endl;
  }
  functionstream << std::string(80, '-');
  log(ERROR, functionstream.str().c_str());
  PyErr_Print();
}

//*******************************************************************|************************************************************//
// initialize the python objects contained in the python instance
//*******************************************************************|************************************************************//
void PythonInstance::init_()
{
  
  PyObject *pGlobals;
  pGlobals = (*GlobalPythonInstance::instance()).globals();
  pLocals_ = PyDict_New();                                           // set up the locals dictionary
  
  pCode_ = PyRun_String((char*)function_.c_str(), Py_file_input,     // run the function string 
                                              pGlobals, pLocals_);
  
  if (PyErr_Occurred()){                                             // check for errors in getting the function
    print_error();
    tf_err("In PythonInstance::init_ evaluating pCode_.", "Python error occurred.");
  }

  pFunc_ = PyDict_GetItemString(pLocals_, "val");                    // get the val function from the function string

  std::stringstream pythonbuffer;                                    // set up a simple python command to check how many 
  pythonbuffer << "import inspect" << std::endl                      // arguments the python function val has
               << "_nargs = len(inspect.getargspec(val).args)" 
               << std::endl;
  PyObject* tmppCode = PyRun_String(pythonbuffer.str().c_str(),      // run the python commands
                                Py_file_input, pGlobals, pLocals_); 
  PyObject* pnArgs = PyDict_GetItemString(pLocals_, "_nargs");       // retrieve the result, _nargs
  nargs_ = PyLong_AsLong(pnArgs);                                     // recast it as an integer
  
  if (PyErr_Occurred()){                                             // check for errors in getting the function
    print_error();
    tf_err("In PythonInstance::init_ evaluating nargs_.", "Python error occurred.");
  }
  
}

//*******************************************************************|************************************************************//
// finalize the python instance
//*******************************************************************|************************************************************//
void PythonInstance::clean_()
{
  Py_DECREF(pLocals_);                                               // decrease the reference count on the locals
  Py_DECREF(pCode_);                                                 // and the code
}

