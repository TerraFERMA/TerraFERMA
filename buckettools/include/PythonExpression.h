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


#ifndef __PYTHON_EXPRESSION_H
#define __PYTHON_EXPRESSION_H

#include <dolfin.h>
#include "Python.h"
#include "PythonInstance.h"
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // PythonExpression class:
  //
  // The PythonExpression class describes a derived dolfin Expression class that overloads
  // the eval function using python
  //*****************************************************************|************************************************************//
  class PythonExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    PythonExpression(const std::string &function);                   // specific constructor (scalar)
    
    PythonExpression(const uint &dim, const std::string &function);  // specific constructor (vector)
    
    PythonExpression(const uint &dim0, const uint &dim1,             // specific constructor (tensor)
                     const std::string &function);
    
    PythonExpression(const std::vector<std::size_t> &value_shape,    // specific constructor (alternate tensor)
                     const std::string &function);
    
    PythonExpression(const std::string &function, 
                                            const double_ptr time);  // specific constructor (scalar)
    
    PythonExpression(const uint &dim, const std::string &function, 
                                            const double_ptr time);  // specific constructor (vector)
    
    PythonExpression(const uint &dim0, const uint &dim1,             // specific constructor (tensor)
                     const std::string &function, const double_ptr time);
    
    PythonExpression(const std::vector<std::size_t> &value_shape,           // specific constructor (alternate tensor)
                     const std::string &function, const double_ptr time);
    
    virtual ~PythonExpression();                                     // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate the expression at a given point
              const dolfin::Array<double>& x) const;                 // (but no cell information?)
    

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//
    
    const bool time_dependent() const;                               // return if this expression is time dependent or not

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//
    
    const PythonInstance pyinst_;                                    // a python instance (wrapping useful python information)

    double_ptr time_;                                                // the time this function is to be evaluated at

  };

}
#endif
