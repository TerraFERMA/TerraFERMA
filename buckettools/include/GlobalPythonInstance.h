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


#ifndef __GLOBAL_PYTHON_INSTANCE_H
#define __GLOBAL_PYTHON_INSTANCE_H

#include "Python.h"
#include <string>

namespace buckettools
{
  //*****************************************************************|************************************************************//
  // GlobalPythonInstance class:
  //
  // The GlobalPythonInstance class wraps global python dictionary so it can be shared between PythonInstances
  //*****************************************************************|************************************************************//
  class GlobalPythonInstance
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    static GlobalPythonInstance* instance();                         // entry point - no public constructor in a singleton class

    void run(const std::string &function);

    PyObject* main()
    { return pMain_; }

    PyObject* globals()
    { return pGlobals_; }

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    GlobalPythonInstance();                                          // specific constructor

    ~GlobalPythonInstance();                                         // default destructor
    

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    PyObject *pMain_, *pGlobals_;                                    // python objects used to run the function (and cacheable between calls)

    static GlobalPythonInstance *instance_;

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
