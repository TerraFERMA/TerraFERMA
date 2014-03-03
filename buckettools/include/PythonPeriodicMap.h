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


#ifndef __PYTHON_PERIODIC_MAP_H
#define __PYTHON_PERIODIC_MAP_H

#include "Python.h"
#include <dolfin.h>
#include "PythonInstance.h"
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // PythonPeriodicMap class:
  //
  // The PythonPeriodicMap class describes a derived dolfin SubDomain class that overloads
  // the map function using python
  //*****************************************************************|************************************************************//
  class PythonPeriodicMap : public dolfin::SubDomain
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    PythonPeriodicMap(const std::string &function);                  // specific constructor
    
    virtual ~PythonPeriodicMap();                                    // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void map(const dolfin::Array<double>& x, dolfin::Array<double>& y) const;        // map slave position x to master position y
    

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//
    
    const PythonInstance pyinst_;                                    // a python instance (wrapping useful python information)

  };

  typedef std::shared_ptr< PythonPeriodicMap > PythonPeriodicMap_ptr;// define a (boost shared) pointer type for the python
                                                                     // periodic map

}
#endif
