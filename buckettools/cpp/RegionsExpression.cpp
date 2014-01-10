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


#include "RegionsExpression.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(std::map< std::size_t, 
                                     Expression_ptr > expressions,
                                     MeshFunction_size_t_ptr cell_ids) : 
                                            dolfin::Expression(), 
                                            expressions_(expressions),
                                            cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const uint &dim, 
                                     std::map< std::size_t, Expression_ptr > 
                                               expressions,
                                     MeshFunction_size_t_ptr cell_ids) : 
                                     dolfin::Expression(dim), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const uint &dim0, const uint &dim1, 
                                     std::map< std::size_t, Expression_ptr > 
                                               expressions,
                                     MeshFunction_size_t_ptr cell_ids) : 
                                     dolfin::Expression(dim0, dim1), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (alternate tensor)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const std::vector<std::size_t>
                                                        &value_shape, 
                                     std::map< std::size_t, Expression_ptr > 
                                               expressions,
                                     MeshFunction_size_t_ptr cell_ids) : 
                                     dolfin::Expression(value_shape), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
RegionsExpression::~RegionsExpression()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// overload dolfin expression eval to loop over a map (from component) of expressions and evaluate them all to a single expression 
// value array
//*******************************************************************|************************************************************//
void RegionsExpression::eval(dolfin::Array<double>& values, 
                                      const dolfin::Array<double>& x, 
                                      const ufc::cell &cell) const
{
  uint id = (*cell_ids_)[cell.index];
  size_t_Expression_const_it e_it = expressions_.find(id);
  if (e_it == expressions_.end())
  {
    dolfin::error("Unknown region id %d in RegionsExpression eval", id);
  }
  else
  {
    (*(*e_it).second).eval(values, x, cell);
  }
}

