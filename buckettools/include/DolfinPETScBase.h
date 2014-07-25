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


#ifndef __DOLFINPETSC_BASE_H
#define __DOLFINPETSC_BASE_H

#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Various tools for converting between dolfin and petsc objects
  //*****************************************************************|************************************************************//

  std::vector<uint> functionspace_dofs(const FunctionSpace_ptr functionspace,
                                       const int &component);

  std::vector<uint> functionspace_dofs_values(const FunctionSpace_ptr functionspace,
                                              MeshFunction_size_t_ptr cellidmeshfunction=NULL,
                                              MeshFunction_size_t_ptr facetidmeshfunction=NULL,
                                              const std::vector<int>* components=NULL,
                                              const std::vector<int>* region_ids=NULL,
                                              const std::vector<int>* boundary_ids=NULL,
                                              PETScVector_ptr values=NULL, 
                                              const dolfin::Expression* value_exp=NULL, const double *value_const=NULL,
                                              uint depth=0, uint exp_index=0);

  boost::unordered_set<uint> cell_dofs_values(const FunctionSpace_ptr functionspace,
                                              MeshFunction_size_t_ptr cellidmeshfunction=NULL,
                                              const std::vector<int>* region_ids=NULL,
                                              PETScVector_ptr values=NULL, 
                                              const dolfin::Expression* value_exp=NULL, const double* value_const=NULL,
                                              const uint &exp_index=0);

  boost::unordered_set<uint> facet_dofs_values(const FunctionSpace_ptr functionspace,
                                               MeshFunction_size_t_ptr facetidmeshfunction=NULL,
                                               const std::vector<int>* boundary_ids=NULL,
                                               PETScVector_ptr values=NULL, 
                                               const dolfin::Expression* value_exp=NULL, const double* value_const=NULL,
                                               const uint &exp_index=0);

  void restrict_indices(std::vector<uint> &indices, 
                        const FunctionSpace_ptr functionspace,
                        const std::vector<uint>* parent_indices=NULL, 
                        const std::vector<uint>* sibling_indices=NULL);

  void restrict_values(PETScVector_ptr values, 
                       PETScVector_ptr tmp_values,
                       const std::vector<uint> &indices);

  IS convert_vector_to_is(const MPI_Comm &comm,
                          const std::vector<uint> &indices,
                          const uint &parent_offset=0, 
                          const std::vector<uint>* parent_indices=NULL);

}

#endif
