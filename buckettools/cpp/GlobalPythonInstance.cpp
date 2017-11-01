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
#include "GlobalPythonInstance.h"
#include "Logger.h"
#include <string>
#include <dolfin.h>

using namespace buckettools;

GlobalPythonInstance* GlobalPythonInstance::instance_ = NULL;        // initialize the global static class variables

//*******************************************************************|************************************************************//
// return the global python instance
//*******************************************************************|************************************************************//
GlobalPythonInstance* GlobalPythonInstance::instance()
{
  if(!instance_)
  {
    instance_ = new GlobalPythonInstance();
  }
  return instance_;
}

//*******************************************************************|************************************************************//
// run a python function from a string
//*******************************************************************|************************************************************//
void GlobalPythonInstance::run(const std::string &function)
{
  PyObject *pCode;
  pCode = PyRun_String((char*)function.c_str(), Py_file_input,     // run the function string 
                                              pGlobals_, pGlobals_);
  
  if (PyErr_Occurred()){                                             // check for errors in getting the function
    log(ERROR, "Global python computation raised an exception.");
    unsigned int lineno = 0;
    std::istringstream functionss(function);
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
    tf_err("In GlobalPythonInstance::run evaluating pCode.", "Python error occurred.");
  }

}

//*******************************************************************|************************************************************//
// specific constructor (private)
//*******************************************************************|************************************************************//
GlobalPythonInstance::GlobalPythonInstance()
{
  init_();                                                           // initialize
}

//*******************************************************************|************************************************************//
// default destructor (private)
//*******************************************************************|************************************************************//
GlobalPythonInstance::~GlobalPythonInstance()
{
  clean_();                                                          // finalize
}

//*******************************************************************|************************************************************//
// initialize the python objects contained in the python instance
//*******************************************************************|************************************************************//
void GlobalPythonInstance::init_()
{
  if(!Py_IsInitialized())
  {
    Py_Initialize();                                                 // if python isn't initialized, do it now
  }
  
  pMain_ = PyImport_AddModule("__main__");                           // set up a main module
  pGlobals_ = PyModule_GetDict(pMain_);                              // set up the globals dictionary
  
}

//*******************************************************************|************************************************************//
// finalize the python instance
//*******************************************************************|************************************************************//
void GlobalPythonInstance::clean_()
{
  if(Py_IsInitialized())                                             // if python is initialized then collect (safe for multiple
  {                                                                  // python instances?)
    PyGC_Collect();
  }
}

