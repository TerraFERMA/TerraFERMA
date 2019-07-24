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


#include <dolfin.h>
#include "Logger.h"
#include "DolfinPETScBase.h"
#include "BucketPETScBase.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// return a vector of dofs from the given functionspace for a field
//*******************************************************************|************************************************************//
std::vector<uint> buckettools::functionspace_dofs(const FunctionSpace_ptr functionspace,
                                                  const int &component)
{
  std::vector<int>* components;
  components = new std::vector<int>;
  (*components).push_back(component);
  std::vector<uint> indices = functionspace_dofs_values(functionspace, NULL, NULL,
                                                        components, NULL, NULL,
                                                        NULL, NULL, NULL);
  delete components;
  components = NULL;
  return indices;
}

//*******************************************************************|************************************************************//
// return a vector of dofs from the given functionspace for a field
//*******************************************************************|************************************************************//
std::vector<uint> buckettools::functionspace_dofs_values(const FunctionSpace_ptr functionspace,
                                                         MeshFunction_size_t_ptr cellidmeshfunction,
                                                         MeshFunction_size_t_ptr facetidmeshfunction,
                                                         const std::vector<int>* components,
                                                         const std::vector<int>* region_ids,
                                                         const std::vector<int>* boundary_ids,
                                                         PETScVector_ptr values, 
                                                         const dolfin::Expression* value_exp, const double *value_const,
                                                         uint depth, uint exp_index)
{
  std::vector<uint> dofs;
  boost::unordered_set<uint> dof_set;

  const uint num_sub_elements = (*(*functionspace).element()).num_sub_elements();
  if (num_sub_elements>0)
  {
    assert(depth==0);                                                // component logic below only makes sense if we only end up
                                                                     // here on the first iteration
    depth++;

    for (uint i = 0; i < num_sub_elements; i++)
    {
      if (components)
      {
        std::vector<int>::const_iterator comp = std::find((*components).begin(), 
                                                          (*components).end(), 
                                                          i);
        if (comp==(*components).end())
        {
          continue;                                                  // component not requested so continue
        }
        exp_index = comp - (*components).begin();                    // work out the index into the expression for this component
      }
      else
      {
        exp_index = i;
      }

      std::vector<uint> tmp_dofs;
      tmp_dofs = functionspace_dofs_values((*functionspace)[i], 
                                           cellidmeshfunction, facetidmeshfunction,
                                           components, region_ids, boundary_ids,
                                           values, value_exp, value_const, 
                                           depth, exp_index);
      dof_set.insert(tmp_dofs.begin(), tmp_dofs.end());
    }

    dofs.insert(dofs.end(), dof_set.begin(), dof_set.end());
    return dofs;
  }
  
  assert(num_sub_elements==0);

  if (boundary_ids)                                                  // do we have boundary id restrictions
  {                                                                  // yes, then get the dofs over these boundaries
    if (region_ids)
    {
      dof_set = cell_dofs_values(functionspace, cellidmeshfunction,  // if we have boundary_ids then we're only interested
                                  region_ids,                        // in cell dofs if we have region_ids specified too
                                  values, value_exp, value_const, 
                                  exp_index);
    }                                                            
    boost::unordered_set<uint> f_dof_set;
    f_dof_set = facet_dofs_values(functionspace, facetidmeshfunction, 
                                  boundary_ids,
                                  values, value_exp, value_const, 
                                  exp_index);
    dof_set.insert(f_dof_set.begin(), f_dof_set.end());
  }
  else                                                               // no boundary_ids specified so let's hope we have some
  {                                                                  // cells to fill the goody bag with
    dof_set = cell_dofs_values(functionspace, cellidmeshfunction, 
                                region_ids,
                                values, value_exp, value_const, 
                                exp_index);
  }

  dofs.insert(dofs.end(), dof_set.begin(), dof_set.end());
  if (values)
  {
    (*values).apply("insert");
  }
  return dofs;

}

//*******************************************************************|************************************************************//
// return a set of dofs from the given functionspace possibly for a subset of the region ids as specified
// FIXME: once mesh domain information is used cellidmeshfunction should be taken directly from the mesh
//*******************************************************************|************************************************************//
boost::unordered_set<uint> buckettools::cell_dofs_values(const FunctionSpace_ptr functionspace,
                                                         MeshFunction_size_t_ptr cellidmeshfunction,
                                                         const std::vector<int>* region_ids,
                                                         PETScVector_ptr values, 
                                                         const dolfin::Expression* value_exp, const double* value_const,
                                                         const uint &exp_index)
{
  boost::unordered_set<uint> dof_set;

  std::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();
  const_Mesh_ptr mesh = (*functionspace).mesh();

  assert((*functionspace).element());
  std::shared_ptr<const dolfin::FiniteElement> element = (*functionspace).element();

  const uint gdim = (*mesh).geometry().dim();                        // set up data for expression evaluation
  boost::multi_array<double, 2> coordinates(boost::extents[(*dofmap).max_cell_dimension()][gdim]);
  std::vector<double> dof_coordinates;
  dolfin::Array<double> x(gdim);

  uint value_size = 1;
  if (value_exp)
  {
    for (uint i = 0; i < (*value_exp).value_rank(); i++)
    {
      value_size *= (*value_exp).value_dimension(i);
    }
    assert(!value_const);
    assert(values);
  }
  dolfin::Array<double> values_array(value_size);

  if (values)
  {
    assert(value_const||value_exp);
  }

  if (value_const)
  {
    assert(!value_exp);
    assert(values);
  }

  for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)       // loop over the cells in the mesh
  {
    if (region_ids)
    {
      int cellid = (*cellidmeshfunction)[(*cell).index()];           // get the cell region id from the mesh function

      std::vector<int>::const_iterator id = std::find((*region_ids).begin(), 
                                             (*region_ids).end(), cellid);
      if (id == (*region_ids).end())
      {
        continue;
      }
    }

    Eigen::Map<const Eigen::Array<dolfin::la_index, Eigen::Dynamic, 1>> tmp_dof_vec = (*dofmap).cell_dofs((*cell).index());
    std::vector<std::size_t> dof_vec;
    for (Eigen::Index i = 0; i < tmp_dof_vec.size(); i++)
    {
      dof_vec.push_back((*dofmap).local_to_global_index(tmp_dof_vec[i]));
    }

    if(value_exp)
    {
      (*cell).get_coordinate_dofs(dof_coordinates);
      (*element).tabulate_dof_coordinates(coordinates, dof_coordinates, *cell);
    }

    for (uint i = 0; i < dof_vec.size(); i++)                        // loop over the cell dof
    {
      dof_set.insert((uint) dof_vec[i]);                             // and insert each one into the unordered set
      if (values)
      {
        if(value_exp)
        {
          for (uint j = 0; j < gdim; j++)
          {
            x[j] = coordinates[i][j];
          }
          (*value_exp).eval(values_array, x);                        // evaluate te expression
          (*values).setitem(dof_vec[i], values_array[exp_index]);    // and set the values to that
        }
        else
        {
          (*values).setitem(dof_vec[i], *value_const);               // and insert each one into the unordered map
                                                                     // assuming a constant
        }
      }
    }
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// return a set of dofs from the given dofmap for the boundary ids specified
// FIXME: once mesh domain information is used facetidmeshfunction should be taken directly from the mesh
//*******************************************************************|************************************************************//
boost::unordered_set<uint> buckettools::facet_dofs_values(const FunctionSpace_ptr functionspace,
                                                          MeshFunction_size_t_ptr facetidmeshfunction,
                                                          const std::vector<int>* boundary_ids,
                                                          PETScVector_ptr values, 
                                                          const dolfin::Expression* value_exp, const double* value_const,
                                                          const uint &exp_index)
{
  boost::unordered_set<uint> dof_set;                                // set up an unordered set of dof

  assert(boundary_ids);

  std::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();
  const_Mesh_ptr mesh = (*functionspace).mesh();

  assert((*functionspace).element());
  std::shared_ptr<const dolfin::FiniteElement> element = (*functionspace).element();

  const uint gdim = (*mesh).geometry().dim();
  boost::multi_array<double, 2> coordinates(boost::extents[(*dofmap).max_cell_dimension()][gdim]);
  std::vector<double> dof_coordinates;
  dolfin::Array<double> x(gdim);

  uint value_size = 1;
  if (value_exp)
  {
    for (uint i = 0; i < (*value_exp).value_rank(); i++)
    {
      value_size *= (*value_exp).value_dimension(i);
    }
    assert(!value_const);
    assert(values);
  }
  dolfin::Array<double> values_array(value_size);

  if (values)
  {
    assert(value_const||value_exp);
  }

  if (value_const)
  {
    assert(!value_exp);
    assert(values);
  }

  for (dolfin::FacetIterator facet(*mesh); !facet.end(); ++facet)   // loop over the facets in the mesh
  {
    int facetid = (*facetidmeshfunction)[(*facet).index()];          // get the facet region id from the mesh function

    for (std::vector<int>::const_iterator id =                       // loop over the region ids that have been requested
                                (*boundary_ids).begin(); 
                                id != (*boundary_ids).end(); id++)
    {
      if(facetid==*id)                                               // check if this facet should be included
      {                                                              // yes...

        const dolfin::Cell cell(*mesh,                               // get cell to which facet belongs
               (*facet).entities((*mesh).topology().dim())[0]);      // (there may be two, but pick first)

        const std::size_t facet_number = cell.index(*facet);         // get the local index of the facet w.r.t. the cell

        Eigen::Map<const Eigen::Array<dolfin::la_index, Eigen::Dynamic, 1>> tmp_cell_dof_vec = (*dofmap).cell_dofs(cell.index());
        std::vector<std::size_t> cell_dof_vec;
        for (Eigen::Index i = 0; i < tmp_cell_dof_vec.size(); i++)
        {
          cell_dof_vec.push_back((*dofmap).local_to_global_index(tmp_cell_dof_vec[i]));
        }

        std::vector<std::size_t> facet_dof_vec((*dofmap).num_facet_dofs(), 0);
        (*dofmap).tabulate_facet_dofs(facet_dof_vec, facet_number);

        if (value_exp)
        {
          cell.get_coordinate_dofs(dof_coordinates);
          (*element).tabulate_dof_coordinates(coordinates, dof_coordinates, cell);
        }

        for (uint i = 0; i < facet_dof_vec.size(); i++)              // loop over facet dof
        {
          dof_set.insert((uint) cell_dof_vec[facet_dof_vec[i]]);     // and insert each one into the unordered set
          if (values)
          {
            if(value_exp)
            {
              for (uint j = 0; j < gdim; j++)
              {
                x[j] = coordinates[i][j];
              }
              (*value_exp).eval(values_array, x);                    // evaluate the values expression
              (*values).setitem((uint) cell_dof_vec[facet_dof_vec[i]], 
                                       values_array[exp_index]);
            }
            else
            {
              (*values).setitem((uint) cell_dof_vec[facet_dof_vec[i]], 
                                       *value_const);                // and insert each one into the unordered map
            }                                                        // assuming a constant
          }
        }                                                         
      }
    }
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// restrict a vector of indices by its parent (intersection) or sibling (complement), also by parallel ownership
//*******************************************************************|************************************************************//
void buckettools::restrict_indices(std::vector<uint> &indices, 
                                   const FunctionSpace_ptr functionspace,
                                   const std::vector<uint>* parent_indices, 
                                   const std::vector<uint>* sibling_indices)
{

  std::vector<uint> tmp_indices = indices;
  indices.clear();

  std::pair<uint, uint> ownership_range =                            // the parallel ownership range of the system functionspace
          (*(*functionspace).dofmap()).ownership_range();

  for (std::vector<uint>::const_iterator                             // loop over the dof in the set
                        dof_it = tmp_indices.begin(); 
                        dof_it != tmp_indices.end(); 
                        dof_it++)
  {                                                                  // and insert them into the indices vector
    if ((*dof_it >= ownership_range.first) &&                        // but first check that this process owns them
                          (*dof_it < ownership_range.second))        // (in parallel)
    {
      indices.push_back(*dof_it);
    }
  }

  std::sort(indices.begin(), indices.end());                         // sort the vector of indices

  if(sibling_indices)                                                // we have been passed a list of sibling indices...
  {                                                                  // we wish to remove from the indices any indices that
                                                                     // also occur in the sibling indices
    tmp_indices.clear();

    std::vector<uint>::const_iterator c_it = indices.begin(),
                                      s_it = (*sibling_indices).begin();
    bool overlap = false;
    while ((c_it != indices.end()) && (s_it != (*sibling_indices).end()))
    {
      if (*c_it < *s_it)
      {
        tmp_indices.push_back(*c_it);
        c_it++;
      }
      else if (*s_it < *c_it)
      {
        s_it++;
      }
      else                                                           // indices[c_ind] == (*sibling_indices)[s_ind]
      {
        c_it++;
        s_it++;
        overlap = true;
      }
    }

    while (c_it != indices.end())                                    // insert any remaining indices beyond the siblings
    {
      tmp_indices.push_back(*c_it);
      c_it++;
    }

    if(overlap)
    {                                                                // sibling indices were ignored... give a warning
      log(WARNING, 
                  "WARNING: IS indices overlap with sibling fieldsplit, ignoring overlapping indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }

  if(parent_indices)                                                 // we have been passed a list of parent indices... 
  {                                                                  // we wish to remove from the indices any indices that do
                                                                     // not occur in the parent indices 
    tmp_indices.clear();

    std::vector<uint>::const_iterator c_it = indices.begin(),
                                      p_it = (*parent_indices).begin();
    bool extra = false;
    while ((c_it != indices.end()) && (p_it != (*parent_indices).end()))
    {
      if (*c_it < *p_it)
      {
        c_it++;                                                      // dropping
        extra = true;
      }
      else if (*p_it < *c_it)
      {
        p_it++;
      }
      else                                                           // *c_it == *p_it
      {
        tmp_indices.push_back(*c_it);
        c_it++;
        p_it++;
      }
    }

    if (c_it != indices.end())
    {
      extra = true;                                                  // dropping trailing
    }

    if(extra)
    {                                                                // indices were ignored... give a warning
      log(WARNING, 
                  "WARNING: IS indices not a subset of parent fieldsplit, ignoring extra indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }
}

void buckettools::restrict_values(PETScVector_ptr values, 
                                  PETScVector_ptr tmp_values,
                                  const std::vector<uint> &indices)
{
  PetscErrorCode perr;                                               // petsc error code

  assert(values);
  assert(tmp_values);

  IS is = convert_vector_to_is((*values).mpi_comm(),                 // this IS indexes from the system vector into whatever the
                               indices);                             // null space vector is (because we haven't supplied
                                                                     // parent_indices for the conversion)
  VecScatter scatter;
  perr = VecScatterCreate((*tmp_values).vec(), is, 
                          (*values).vec(), is, 
                          &scatter);
  petsc_err(perr);
  perr = VecScatterBegin(scatter, 
                         (*tmp_values).vec(), (*values).vec(), 
                         INSERT_VALUES, SCATTER_FORWARD);
  petsc_err(perr);
  perr = VecScatterEnd(scatter,
                       (*tmp_values).vec(), (*values).vec(),
                       INSERT_VALUES, SCATTER_FORWARD);
  petsc_err(perr);

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1            // necessary or taken care of when object leaves scope?
  perr = VecScatterDestroy(&scatter); petsc_err(perr);      
  perr = ISDestroy(&is); petsc_err(perr);
  #else
  perr = VecScatterDestroy(scatter); petsc_err(perr);
  perr = ISDestroy(is); petsc_err(perr);
  #endif

}

//*******************************************************************|************************************************************//
// Fill a petsc IS object based on the supplied indices.
//*******************************************************************|************************************************************//
IS buckettools::convert_vector_to_is(const MPI_Comm &comm,
                                     const std::vector<uint> &indices,
                                     const uint &parent_offset, 
                                     const std::vector<uint>* parent_indices)
{

  PetscErrorCode perr;                                               // petsc error code

  std::vector<uint>::const_iterator pos =  std::adjacent_find(indices.begin(), 
                                                              indices.end(),  
                                                              std::greater<uint>());
  assert(pos==indices.end());                                        // test that the incoming indices are sorted

  PetscInt n=indices.size();                                         // setup a simpler structure for petsc
  if (dolfin::MPI::size(comm) == 1)
  {
    assert(n>0);                                                     // not necessarily true on more than one process
  }
  PetscInt *pindices;
  PetscMalloc(n*sizeof(PetscInt), &pindices);

  uint ind = 0;
  if(parent_indices)
  {                                                                  // we have been passed a list of parent indices... 
                                                                     // our child indices must be a  subset of this list and indexed
                                                                     // into it so let's do that now while we convert structures...
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        it = indices.begin(); 
                                        it != indices.end(); 
                                        it++)
    {
      while ((*parent_indices)[p_ind] != *it)                      // indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and throw an error
          tf_err("IS indicies are not a subset of a parent fieldsplit.", "p_ind = %d, p_size = %d", p_ind, p_size);
        }
      }
      pindices[ind] = p_ind + parent_offset;                         // found the child index in the parent_indices so copy it into
                                                                     // the PetscInt array
      ind++;                                                         // increment the array index
      p_ind++;                                                       // indices shouldn't be repeated so increment the parent too
    } 
    assert(ind==n);                                                  // these should be equal
  }
  else
  {
    for (std::vector<uint>::const_iterator                           // loop over the indices
                                      ind_it = indices.begin(); 
                                      ind_it != indices.end(); 
                                      ind_it++)
    {
      pindices[ind] = *ind_it;                                       // insert them into the PetscInt array
      ind++;                                                         // increment the array index
    }
    assert(ind==n);                                                  // these should be equal
  }

  IS is;
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  perr = ISCreateGeneral(comm, n, pindices, PETSC_OWN_POINTER, &is); // create the general index set based on the indices
  #else
  perr = ISCreateGeneral(comm, n, pindices, &is);                    // create the general index set based on the indices
  #endif

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  #else
  PetscFree(pindices);                                               // free the PetscInt array of indices
  #endif

  return is;
}



