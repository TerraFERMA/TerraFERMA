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


#include "InitialConditionExpression.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
InitialConditionExpression::InitialConditionExpression(std::map< std::size_t, 
                                     Expression_ptr > expressions) : 
                                            dolfin::Expression(), 
                                            expressions_(expressions)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
InitialConditionExpression::InitialConditionExpression(const uint &dim, 
                      std::map< std::size_t, Expression_ptr > expressions) : 
                      dolfin::Expression(dim), 
                      expressions_(expressions)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
InitialConditionExpression::~InitialConditionExpression()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// overload dolfin expression eval to loop over a map (from component) of expressions and evaluate them all to a single expression 
// value array
//*******************************************************************|************************************************************//
void InitialConditionExpression::eval(dolfin::Array<double>& values, 
                                      const dolfin::Array<double>& x, 
                                      const ufc::cell &cell) const
{
  for (uint i = 0; i < values.size(); i++)                           // zero the full value array
  {
    values[i] = 0.0;
  }
  for (size_t_Expression_const_it expr = expressions_.begin();         // loop over the expressions in the map
                    expr != expressions_.end(); expr++)
  {
    dolfin::Array<double> tmp_values((*(*expr).second).value_size());// set up a temporary value array
    (*(*expr).second).eval(tmp_values, x, cell);                     // evaluate the component expression
    for (uint i = 0; i< tmp_values.size(); i++)                      // loop over the temporary value array
    {
      values[(*expr).first+i] = tmp_values[i];                       // insert component value into the full value array
    }
  }
  
}
