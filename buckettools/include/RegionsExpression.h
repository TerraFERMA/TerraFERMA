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


#ifndef __REGIONS_EXPRESSION_H
#define __REGIONS_EXPRESSION_H

#include <dolfin.h>
#include "BoostTypes.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // RegionsExpression class:
  //
  // This class provides a method of overloading a dolfin Expression eval function by looping over cell ids and returning the results
  // of individual expressions in each of those regions.
  //*****************************************************************|************************************************************//
  class RegionsExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    RegionsExpression(                                               // specific constructor (scalar)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const uint &dim,                               // specific constructor (vector)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const uint &dim0, const uint &dim1,            // specific constructor (tensor)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const std::vector<size_t> &value_shape,          // specific constructor (alternate tensor)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    virtual ~RegionsExpression();                           // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate this expression at a point in a given cell
              const dolfin::Array<double>& x, 
              const ufc::cell &cell) const;
    
    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::map< std::size_t, Expression_ptr> expressions() const        // return the map of expressions (const version)
    { return expressions_; }
    
    std::map< std::size_t, Expression_ptr> expressions()                    // return the map of expressions (non-const version)
    { return expressions_; }
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::size_t, Expression_ptr > expressions_;                   // map from region id to expression for that region
  
    MeshFunction_size_t_ptr cell_ids_;                                 // a (boost shared) pointer to a (uint) mesh function holding the cell ids

  };

}
#endif
