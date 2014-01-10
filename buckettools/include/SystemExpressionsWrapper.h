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

#ifndef __SYSTEMEXPRESSIONS_WRAPPER_H
#define __SYSTEMEXPRESSIONS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  Expression_ptr cpp_fetch_expression(const std::string &systemname, // return a (boost shared) pointer to an expression 
                          const std::string &functionname,           // given a system, function & expression name and a bucket
                          const std::string &expressiontype,         // and an expression type (initial_condition, boudary_condition,
                          const std::string &expressionname,         // value)
                          const uint &size, 
                          const std::vector<int> &shape, 
                          const Bucket *bucket,  
                          const SystemBucket *system,
                          const double_ptr time);

  void cpp_init_expression(Expression_ptr expression,                // initialize the given (boost shared) pointer to an expression
                          const std::string &systemname,             // given a system, function & expression name & an expression type
                          const std::string &functionname,
                          const std::string &expressiontype,
                          const std::string &expressionname);

}

#endif
