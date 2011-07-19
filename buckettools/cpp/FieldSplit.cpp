// Copyright (C) 2011 Marc Spiegelman 
// Licensed under the GNU LGPL Version 3.1.
//
// First added:  16 Jun 2011 11:45:10
// Last changed:  16 Jun 2011 11:45:15
//
// set of utility routines for converting and combining Dolfin Degrees of freedom
// into PETSc Index Sets (IS)
// 
// 

#include "petsc.h"
#include "petscis.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include "FieldSplit.h"


using namespace dolfin;

//----------------------------------------------------------------------------
// Constructor:
//----------------------------------------------------------------------------
FieldSplit::FieldSplit(const FunctionSpace& V) : V(V), _num_fields(0) 
{
  _field_dofs = get_field_dofs();
  _num_fields = _field_dofs.size();
  _ghost_indices = get_ghost_indices();

}
//----------------------------------------------------------------------------
// Destructor:
//----------------------------------------------------------------------------
FieldSplit::~FieldSplit()
{
  // do nothing
};

//----------------------------------------------------------------------------
// return number of fields field dofs:
//----------------------------------------------------------------------------
uint FieldSplit::num_fields() const
{
  return _num_fields;
};

//----------------------------------------------------------------------------
// return full set of field dofs:
//----------------------------------------------------------------------------
std::vector<boost::unordered_set<dolfin::uint> > FieldSplit::field_dofs() const
{
  return _field_dofs;
};

//----------------------------------------------------------------------------
// return  field dofs for field i;
//----------------------------------------------------------------------------
boost::unordered_set<dolfin::uint> FieldSplit::field_dofs(const uint i) const
{
  return _field_dofs[i];
};

//----------------------------------------------------------------------------
// return  total ghost indices dofs for field i;
//----------------------------------------------------------------------------
std::vector<uint> FieldSplit::ghost_indices() const
{
  return _ghost_indices;
};

//-----------------------------------------------------------------------------
//combine and sort dofmaps for specific splits and load a Petsc Index Set
//-----------------------------------------------------------------------------
IS  FieldSplit::getIS(std::vector<dolfin::uint>  &split) 
{
  // output sorted vector
  std::vector<dolfin::uint> isv;
  
  //split iterator
  std::vector<dolfin::uint>::iterator it;
  //dof iterator
  boost::unordered_set<dolfin::uint> dof_i;
  boost::unordered_set<dolfin::uint>::iterator dof_it;

  //loop over split sub_fields and concatenate vector
  for (it = split.begin(); it != split.end(); it++)
  {
    dof_i = _field_dofs[*it];
    for (dof_it = dof_i.begin(); dof_it != dof_i.end(); dof_it++)
    {
      isv.push_back(*dof_it);
    }
  }

  //sort vector
  std::sort(isv.begin(),isv.end());
  
 
  //allocate maximum storage for IS indices
  PetscInt n=isv.size();
  PetscInt *indices;
  
  PetscMalloc(n*sizeof(PetscInt),&indices);

  //load indices restricted to ownership range
  std::pair<uint, uint> ownership_range = V.dofmap().ownership_range();

  uint i=0;
  for (it = isv.begin(); it != isv.end(); it++)
  {
    if ((*it >= ownership_range.first) && (*it < ownership_range.second)) 
    {
      indices[i] = *it;
      i++;
    }
  }

  //set count size to i;
  n = i;
  //Create General IS
  IS is;
  ISCreateGeneral(PETSC_COMM_WORLD,n,indices,&is);

  //free indices
  PetscFree(indices);

  return is;

};

//-----------------------------------------------------------------------------
//utility function to print degrees of freedom for each Field
//-----------------------------------------------------------------------------
void FieldSplit::print_dofs() 
{
  
  boost::unordered_set<dolfin::uint>::iterator it;
  
  for (uint i=0; i< _num_fields; i++ )
  {
    cout << "Field[" << i << "] dofs: ";
    for (it = _field_dofs[i].begin(); it != _field_dofs[i].end(); it++ )
    {
      cout << *it << "," ;
    }
    cout << endl;
  }
};

//-----------------------------------------------------------------------------
// get_field_dofs: Function to extract degrees of freedom indices for each sub-field of a Function
//-----------------------------------------------------------------------------
std::vector<boost::unordered_set<dolfin::uint> >  FieldSplit::get_field_dofs()
{
  //allocate vector 
  std::vector<boost::unordered_set<dolfin::uint> > field_dofs;

  //loop over sub-elements allowing two levels of sub-spaces
  //fix me, there must be a cleaner way to do this

  for (dolfin::uint i=0; i< V.element().num_sub_elements(); i++)
  {
    boost::shared_ptr<dolfin::FunctionSpace> Vs, Vss;
    Vs = V[i];
    dolfin::uint n_sub_subs = Vs->element().num_sub_elements();
    if (n_sub_subs == 0)
    {
      field_dofs.push_back(Vs->dofmap().dofs());
    }
    else
    {
      for (dolfin::uint j=0; j<n_sub_subs; j++) 
      {
        Vss = (*Vs)[j];
        field_dofs.push_back(Vss->dofmap().dofs());
      }
    } 
  }
  return field_dofs;
};

 
//-----------------------------------------------------------------------------
// get_ghost_indices: Function to extract degrees of freedom in function space, that are not in ownership range
//-----------------------------------------------------------------------------
std::vector<uint>  FieldSplit::get_ghost_indices()
{
  //allocate vector 
  std::vector<uint> ghost_indices;

  //get all dofs for function space

  //fix me, there must be a cleaner way to do this

  boost::unordered_set<dolfin::uint> dofs = V.dofmap().dofs();
  boost::unordered_set<dolfin::uint>::iterator dof_it;

  //get global size
  const uint N = V.dofmap().global_dimension();

  //get local range
  const std::pair<uint, uint> range = V.dofmap().ownership_range();
  const uint n0 = range.first;
  const uint n1 = range.second;
  const uint local_size = n1 - n0;

  //if distributed
  if (N > local_size)
  {
    for (dof_it = dofs.begin(); dof_it != dofs.end(); dof_it++ )
    {
      if (*dof_it < n0 || *dof_it >= n1 )
      {
        ghost_indices.push_back(*dof_it);
      }
    }
    std::sort(ghost_indices.begin(),ghost_indices.end());
  }

  return ghost_indices;
};
