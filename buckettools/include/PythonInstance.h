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


#ifndef __PYTHON_INSTANCE_H
#define __PYTHON_INSTANCE_H

#include "Python.h"
#include <string>

namespace buckettools
{
  //*****************************************************************|************************************************************//
  // PythonInstance class:
  //
  // The PythonInstance class wraps python functionality so that it is accessible from C++
  // NOTE: This assumes a function name val is available and runs that and only that!
  //*****************************************************************|************************************************************//
  class PythonInstance
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    PythonInstance(const std::string &function);                     // specific constructor (takes a string with the python function)
    
    ~PythonInstance();                                               // default destructor
    
    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string function() const                               // return a constant string containing the function
    { return function_; }
    
    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    PyObject* call(PyObject *pArgs) const;                           // run the function contained in this python instance

    const int number_arguments() const                               // return the number of arguments expected by this pythoninstance
    { return nargs_; }

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    const std::string function_;                                     // the python function string

    PyObject *pMain_, *pGlobals_, *pLocals_, *pCode_, *pFunc_;       // python objects used to run the function (and cacheable between calls)

    int nargs_;                                                      // the number of arguments this python function takes
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_();                                                    // initialize the python instance
    
    //***************************************************************|***********************************************************//
    // Clean up
    //***************************************************************|***********************************************************//

    void clean_();                                                   // clean the python instance
    
  };
}
#endif
