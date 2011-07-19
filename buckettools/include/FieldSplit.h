// Copyright (C) 2011 Marc Spiegelman 
// Licensed under the GNU LGPL Version 3.1.
//
// First added:  16 Jun 2011 11:45:10
// Last changed:  20 Jun 2011 13:34:05
//
// Class for converting Dolfin Degrees of freedom and subfunctions to PETSc IS's for fieldsplit preconditioners
// 
// 

#ifndef __FIELDSPLIT_H
#define __FIELDSPLIT_H

#include "petsc.h"
#include "petscis.h"
#include "BoostTypes.h"
#include <dolfin.h>


namespace buckettools
{

  // Forward declarations
  class FunctionSpace;
  class DofMap;

  // This class allows conversion between Dolfin Degrees of freedom on FunctionSpaces and sub-function spaces to
  // PETSc Index Sets (IS's)  for eventual use in defining FieldSplit block preconditioners from Dolfin Functions

  class FieldSplit 
  {
  public:

    // Create FieldSplit for a given function space
    explicit FieldSplit(const FunctionSpace& V);

    // Destructor
    virtual ~FieldSplit();

    // return number of sub fields
    uint num_fields() const;

    // return field degrees of freedom
    std::vector<boost::unordered_set<dolfin::uint> > field_dofs() const;

    // return field degrees of freedom for field i
    boost::unordered_set<dolfin::uint>  field_dofs(const uint i) const;    

    // return Index Set for specific split
    IS getIS(std::vector<uint>& split);

    // return ownership range
    std::pair<uint,uint> range() { return V.dofmap().ownership_range(); }

    // return ghost indices
    std::vector<uint> ghost_indices() const;


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
    uint _num_fields;

    // vector of field_dofs
    mutable std::vector<boost::unordered_set<dolfin::uint> > _field_dofs;

    //vector of ghost_indices (all dofs not in ownership range)
    std::vector<uint> _ghost_indices;

    // function to extract field dofs
    std::vector<boost::unordered_set<dolfin::uint> >  get_field_dofs();

    // function to extract ghost indices
    std::vector<uint>  get_ghost_indices();  
  };
}

#endif
