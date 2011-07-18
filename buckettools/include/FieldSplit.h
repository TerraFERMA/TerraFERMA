// Copyright (C) 2011 Marc Spiegelman 
// Licensed under the GNU LGPL Version 3.1.
//
// First added:  16 Jun 2011 11:45:10
// Last changed:  20 Jun 2011 13:34:05
//
// Class for converting Dolfin Degrees of freedom and subfunctions to PETSc IS's for fieldsplit preconditioners
// 
// 

#ifndef __DOLFINFIELDSPLIT_H
#define __DOLFINFIELDSPLIT_H

#include "petsc.h"
#include "petscis.h"
/* #include <dolfin/fem/DofMap.h> */
/* #include <dolfin/function/Function.h> */
/* #include <dolfin/function/FunctionSpace.h> */
#include <boost/unordered_set.hpp>
#include <dolfin/common/types.h>
#include <dolfin/fem/DofMap.h>

namespace dolfin
{

  // Forward declarations
  class FunctionSpace;
  class DofMap;

  // This class allows conversion between Dolfin Degrees of freedom on FunctionSpaces and sub-function spaces to
  // PETSc Index Sets (IS's)  for eventual use in defining FieldSplit block preconditioners from Dolfin Functions

  class DolfinFieldSplit 
  {
  public:

    // Create DolfinFieldSplit for a given function space
    explicit DolfinFieldSplit(const FunctionSpace& V);

    // Destructor
    virtual ~DolfinFieldSplit();

    // return number of sub fields
    unsigned int num_fields() const;

    // return field degrees of freedom
    std::vector<boost::unordered_set<dolfin::uint> > field_dofs() const;

    // return field degrees of freedom for field i
    boost::unordered_set<dolfin::uint>  field_dofs(const unsigned int i) const;    

    // return Index Set for specific split
    IS getIS(std::vector<unsigned int>& split);

    // return ownership range
    std::pair<unsigned int,unsigned int> range() { return V.dofmap().ownership_range(); }

    // return ghost indices
    std::vector<unsigned int> ghost_indices() const;


    // set SNES Fieldsplit (where Fieldsplit is a tree of splits
    // void set(SNES snes, FieldSplit fieldsplit);

    //set KSP Fieldsplit
    //void set(KSP ksp, FieldSplit fieldsplit);

    // utility function to print subfunction dofs
    void print_dofs();

    

  private:
    
    // Function space
    const FunctionSpace& V;
    
    //number of sub_fields
    unsigned int _num_fields;

    // vector of field_dofs
    mutable std::vector<boost::unordered_set<dolfin::uint> > _field_dofs;

    //vector of ghost_indices (all dofs not in ownership range)
    std::vector<unsigned int> _ghost_indices;

    // function to extract field dofs
    std::vector<boost::unordered_set<dolfin::uint> >  get_field_dofs();

    // function to extract ghost indices
    std::vector<unsigned int>  get_ghost_indices();  
  };
}

#endif
