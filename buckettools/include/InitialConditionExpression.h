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


#ifndef __INITIAL_CONDITION_EXPRESSION_H
#define __INITIAL_CONDITION_EXPRESSION_H

#include "Bucket.h"
#include <dolfin.h>
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // InitalConditionExpression class:
  //
  // This class provides a method of looping over field initial conditions (as individual expressions)
  // and setting a mixed function space to those initial conditions.
  // It is a derived dolfin expression and overloads much functionality in that base class
  //*****************************************************************|************************************************************//
  class InitialConditionExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    InitialConditionExpression(                                      // specific constructor (scalar)
                      std::map< std::size_t, Expression_ptr > expressions);
    
    InitialConditionExpression(const uint &dim,                      // specific constructor (vector)
                      std::map< std::size_t, Expression_ptr > expressions);
    
    InitialConditionExpression(const std::vector<std::size_t> &value_shape,// specific constructor (tensor)
                      std::map< std::size_t, Expression_ptr > expressions);
    
    
    virtual ~InitialConditionExpression();                           // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate this expression at a point in a given cell
              const dolfin::Array<double>& x, 
              const ufc::cell &cell) const;
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::size_t, Expression_ptr > expressions_;                   // map from component to initial condition expression for a function
  
  };

}
#endif
