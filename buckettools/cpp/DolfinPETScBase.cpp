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
#include "DolfinPETScBase.h"

using namespace buckettools;

std::vector<uint> buckettools::functionspace_dofs(const FunctionSpace_ptr functionspace,
                                                  const int &component)
{
  std::vector<int>* components;
  components = new std::vector<int>;
  (*components).push_back(component);
  std::vector<uint> indices = functionspace_dofs(functionspace, NULL, NULL,
                                                 components, NULL, NULL);
  delete components;
  components = NULL;
  return indices;
}
//*******************************************************************|************************************************************//
// return a vector of dofs from the given functionspace for a field
//*******************************************************************|************************************************************//
std::vector<uint> buckettools::functionspace_dofs(const FunctionSpace_ptr functionspace,
                                                  MeshFunction_size_t_ptr cellidmeshfunction,
                                                  MeshFunction_size_t_ptr facetidmeshfunction,
                                                  const std::vector<int>* components,
                                                  const std::vector<int>* region_ids,
                                                  const std::vector<int>* boundary_ids,
                                                  uint depth)
{
  assert(depth<=1);                                                   // component logic below only makes sense if we've 
                                                                      // only iterated at most once

  std::vector<uint> dofs;
  boost::unordered_set<uint> dof_set;

  const uint num_sub_elements = (*(*functionspace).element()).num_sub_elements();
  if (num_sub_elements>0)
  {
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
      }

      std::vector<uint> tmp_dofs;
      tmp_dofs = functionspace_dofs((*functionspace)[i], 
                                    cellidmeshfunction, facetidmeshfunction,
                                    components, region_ids, boundary_ids,
                                    depth);
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
      dof_set = cell_dof_set(functionspace, cellidmeshfunction, 
                                            region_ids);             // if we have boundary_ids then we're only interested
    }                                                                // in cell dofs if we have region_ids specified too
    boost::unordered_set<uint> f_dof_set;
    f_dof_set = facet_dof_set(functionspace, facetidmeshfunction, 
                                                 boundary_ids);
    dof_set.insert(f_dof_set.begin(), f_dof_set.end());
  }
  else                                                               // no boundary_ids specified so let's hope we have some
  {                                                                  // cells to fill the goody bag with
    dof_set = cell_dof_set(functionspace, cellidmeshfunction, 
                                          region_ids);
  }

  dofs.insert(dofs.end(), dof_set.begin(), dof_set.end());
  return dofs;

}

//*******************************************************************|************************************************************//
// return a set of dofs from the given functionspace possibly for a subset of the region ids as specified
// FIXME: once mesh domain information is used cellidmeshfunction should be taken directly from the mesh
//*******************************************************************|************************************************************//
boost::unordered_set<uint> buckettools::cell_dof_set(const FunctionSpace_ptr functionspace,
                                                     MeshFunction_size_t_ptr cellidmeshfunction,
                                                     const std::vector<int>* region_ids)
{
  boost::unordered_set<uint> dof_set;

  std::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();
  const_Mesh_ptr mesh = (*functionspace).mesh();

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

    std::vector<dolfin::la_index> dof_vec = (*dofmap).cell_dofs((*cell).index());
    for (std::vector<dolfin::la_index>::const_iterator dof_it =      // loop over the cell dof
                                    dof_vec.begin(); 
                                    dof_it < dof_vec.end(); 
                                    dof_it++)
    {
      dof_set.insert((uint) *dof_it);                                // and insert each one into the unordered set
    }                                                                // (i.e. if it hasn't been added already)
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// return a set of dofs from the given dofmap for the boundary ids specified
// FIXME: once mesh domain information is used facetidmeshfunction should be taken directly from the mesh
//*******************************************************************|************************************************************//
boost::unordered_set<uint> buckettools::facet_dof_set(const FunctionSpace_ptr functionspace,
                                                      MeshFunction_size_t_ptr facetidmeshfunction,
                                                      const std::vector<int>* boundary_ids)
{
  boost::unordered_set<uint> dof_set;                                // set up an unordered set of dof

  std::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();
  const_Mesh_ptr mesh = (*functionspace).mesh();

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

        std::vector<dolfin::la_index> cell_dof_vec;
        cell_dof_vec = (*dofmap).cell_dofs(cell.index());            // get the cell dof (potentially for all components)
        
        std::vector<std::size_t> facet_dof_vec((*dofmap).num_facet_dofs(), 0);
        (*dofmap).tabulate_facet_dofs(facet_dof_vec, facet_number);

        for (std::vector<std::size_t>::const_iterator dof_it =       // loop over the cell dof
                                facet_dof_vec.begin(); 
                                dof_it < facet_dof_vec.end(); 
                                dof_it++)
        {
          dof_set.insert((uint) cell_dof_vec[*dof_it]);              // and insert each one into the unordered set
        }                                                            // (i.e. if it hasn't been added already)
      }
    }
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// restrict a vector of indices by its parent (intersection) or sibling (complement), also by parallel ownership
//*******************************************************************|************************************************************//
void buckettools::restrict_is_indices(std::vector<uint> &indices, 
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

    uint c_size = indices.size();
    uint c_ind = 0;
    bool overlap = false;
    for(std::vector<uint>::const_iterator                            // loop over the sibling indices
                                   s_it = (*sibling_indices).begin();
                                   s_it != (*sibling_indices).end();
                                   s_it++)
    {
      while(indices[c_ind] != *s_it)                                 // indices are sorted, so sibling_indices should be too
      {                                                              // search indices until the current sibling index is found
        tmp_indices.push_back(indices[c_ind]);                       // include indices that aren't in the sibling
        c_ind++;
        if (c_ind == c_size)                                         // or we reach the end of the indices...
        {
          break;
        }
      }
      if (c_ind == c_size)                                           // we've reached the end of the indices so nothing more
      {                                                              // to do
        break;
      }
      else                                                           // we haven't reached the end of the indices but found
      {                                                              // a sibling index to ignore... give a warning
        overlap = true;
      }
      c_ind++;                                                       // indices shouldn't be repeated so increment the too
    }
    for (uint i = c_ind; i < c_size; i++)
    {
      tmp_indices.push_back(indices[i]);                             // insert any remaining indices beyond the siblings
    }

    if(overlap)
    {                                                                // sibling indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices overlap with sibling fieldsplit, ignoring overlapping indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }

  if(parent_indices)                                                 // we have been passed a list of parent indices... 
  {                                                                  // we wish to remove from the indices any indices that do
                                                                     // not occur in the parent indices 
    tmp_indices.clear();

    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    uint p_reset = 0;
    bool extra = false;
    for (std::vector<uint>::const_iterator                           // loop over the indices
                                        c_it = indices.begin(); 
                                        c_it != indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and prepare to throw a warning
          extra = true;
          break;
        }
      }
      if (p_ind == p_size)
      {
        p_ind = p_reset;
      }
      else
      {
        tmp_indices.push_back(*c_it);                                // include indices that are in the parent
        p_ind++;                                                     // indices shouldn't be repeated so increment the parent too
        p_reset = p_ind;                                             // this is where the next failed search should continue from
        if (p_ind == p_size)                                         // we've reached the end
        { 
          break;                                            
        }
      }
    } 

    if(extra)
    {                                                                // indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices not a subset of parent fieldsplit, ignoring extra indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }
}

//*******************************************************************|************************************************************//
// Fill a petsc IS object based on the supplied indices.
//*******************************************************************|************************************************************//
IS buckettools::convert_vector_to_is(const MPI_Comm &comm,
                                     const std::vector<uint> &indices,
                                     const uint &offset, 
                                     const std::vector<uint>* parent_indices)
{

  PetscErrorCode perr;                                               // petsc error code

  std::vector<uint>::const_iterator pos =  std::adjacent_find(indices.begin(), 
                                                              indices.end(),  
                                                              std::greater<uint>());
  assert(pos==indices.end());                                        // test that the incoming indices are sorted

  PetscInt n=indices.size();                                         // setup a simpler structure for petsc
  assert(n>0);
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
          dolfin::error("IS indices are not a subset of a parent fieldsplit, shouldn't happen here.");
        }
      }
      pindices[ind] = p_ind + offset;                                // found the child index in the parent_indices so copy it into
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



