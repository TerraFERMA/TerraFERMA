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

#ifndef __SYSTEMSOLVERS_WRAPPER_H
#define __SYSTEMSOLVERS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  FunctionSpace_ptr ufc_fetch_functionspace(const std::string        // return a (boost shared) pointer to a functionspace from a
                                            &systemname,             // system given a (boost shared) pointer to a mesh (defaults to the
                                      Mesh_ptr mesh);                // first solver in the system as they should all be the same)

  FunctionSpace_ptr ufc_fetch_functionspace(const std::string        // return a (boost shared) pointer to a functionspace from a
                                            &systemname,             // system given a mesh and a solver name (i.e. does not default
                                      const std::string &solvername, // as above) 
                                      Mesh_ptr mesh);

  FunctionSpace_ptr ufc_fetch_coefficientspace_from_solver(
                                      const std::string              // return a (boost shared) pointer to a functionspace for a
                                            &systemname,             // coefficient function from a solver given a mesh, a solver name
                                      const std::string &solvername, // and a (base) ufl symbol
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh);

  Form_ptr ufc_fetch_form(const std::string &systemname,             // return a (boost shared) pointer to a form from a solver
                          const std::string &solvername,             // given a functionspace, a solver name, a solver type
                          const std::string &solvertype,             // and a form name
                          const std::string &formname, 
                          FunctionSpace_ptr functionspace);

}

#endif
