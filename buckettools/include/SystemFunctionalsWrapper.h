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

#ifndef __SYSTEMFUNCTIONALS_WRAPPER_H
#define __SYSTEMFUNCTIONALS_WRAPPER_H

#include "PythonPeriodicMap.h"
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(
                                      const std::string              // return a (boost shared) pointer to a functionspace for a 
                                            &systemname,             // coefficient function from a functional given a mesh, a
                                      const std::string              // functional name and a (base) ufl symbol
                                            &funcionalname, 
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh,
                                      PythonPeriodicMap_ptr periodicmap,
                                      MeshFunction_size_t_ptr facetdomains,
                                      const std::vector<std::size_t> &masterids,
                                      const std::vector<std::size_t> &slaveids);

  FunctionSpace_ptr ufc_fetch_coefficientspace_from_constant_functional(
                                      const std::string              // return a (boost shared) pointer to a functionspace for a 
                                            &systemname,             // coefficient function from a constant functional given a mesh, a
                                      const std::string              // function name and a (base) ufl symbol
                                            &coefficientname, 
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh,
                                      PythonPeriodicMap_ptr periodicmap,
                                      MeshFunction_size_t_ptr facetdomains,
                                      const std::vector<std::size_t> &masterids,
                                      const std::vector<std::size_t> &slaveids);

  Form_ptr ufc_fetch_functional(const std::string &systemname,       // return a (boost shared) pointer to a form from a functional
                          const std::string &functionalname,         // given a mesh and a functional name
                          Mesh_ptr mesh);

  Form_ptr ufc_fetch_constant_functional(const std::string &systemname,// return a (boost shared) pointer to a form for a constant
                          const std::string &functionname,           // from a functional given a mesh and a function name
                          Mesh_ptr mesh);

}

#endif
