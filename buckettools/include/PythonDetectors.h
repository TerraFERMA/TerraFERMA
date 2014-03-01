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


#ifndef __PYTHON_DETECTORS_H
#define __PYTHON_DETECTORS_H

#include "PythonInstance.h"
#include "GenericDetectors.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // PythonDetectors class:
  //
  // PythonDetectors is a derived class of GenericDetectors that implements an array of detectors described by a python function
  //*****************************************************************|************************************************************//
  class PythonDetectors : public GenericDetectors
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    PythonDetectors(const uint &meshdim, 
                    const std::string &function, 
                    const std::string &name);
    
    ~PythonDetectors();
    
    //***************************************************************|***********************************************************//
    // Output
    //***************************************************************|***********************************************************//

    const std::string str() const;                                 // return a string that describes the detectors 
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    PythonInstance pyinst_;                                          // an object that describes the base python function
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_();                                                    // initialize the python detectors
    
  };
  
  typedef std::shared_ptr< PythonDetectors > PythonDetectors_ptr;  // define a (boost shared) pointer for this class type
  
}

#endif
