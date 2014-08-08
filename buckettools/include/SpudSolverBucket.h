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


#ifndef __SPUD_SOLVERBUCKET_H
#define __SPUD_SOLVERBUCKET_H

#include "BoostTypes.h"
#include "SolverBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudSolverBucket class:
  //
  // The SpudSolverBucket class is a derived class of the solver that populates the
  // data structures within a solver using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudSolverBucket : public SolverBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudSolverBucket(const std::string &optionpath, SystemBucket* system);  // specific constructor (taking in optionpath and parent system)
    
    ~SpudSolverBucket();                                             // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill in the data in the base solver bucket

    void initialize();                                               // initialize the solvers and tensors

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a constant string containing the optionpath for the solver bucket
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Form data access
    //***************************************************************|***********************************************************//

    void register_form(Form_ptr form, const std::string &name,       // register a form with the given name and optionpath in the solver
                                   std::string optionpath);   // bucket

    const std::string fetch_form_optionpath(const std::string &name) const;// return the optionpath of the named form

    string_it form_optionpaths_begin();                          

    string_const_it form_optionpaths_begin() const;   

    string_it form_optionpaths_end();              

    string_const_it form_optionpaths_end() const;   

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str(int indent=0) const;                       // return an indented string describing the contents of the solver bucket

    const std::string forms_str(const int &indent=0) const;           // return an indented string describing the forms in the solver bucket

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible by this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the optionpath for the solver bucket

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    ordered_map< const std::string, std::string > form_optionpaths_;          // a map from form names to form optionpaths
    
    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill the base data of the solver bucket
 
    void fill_forms_();                                              // fill the form data of the solver bucket

    void fill_subforms_(const std::string &optionpath, 
                        const std::string &prefix="");               // fill the form data of a section of the solver bucket

    void fill_solverforms_(const std::string &optionpath, 
                           const std::string &prefix="");            // fill the form data of a linear solver

    void fill_ksp_(const std::string &optionpath, KSP &ksp,          // fill the information about a child ksp
                   const std::string prefix, 
                   const std::vector<uint>* parent_indices=NULL);

    void fill_pc_(const std::string &optionpath, PC &pc,             // fill the information about a pc
                  const std::string prefix, 
                  const uint &parent_offset, 
                  const std::vector<uint>* parent_indices);

    void fill_indices_values_by_field_(const std::string &optionpath,
                                       std::vector<uint> &child_indices,
                                       PETScVector_ptr values=NULL,
                                       const std::vector<uint>* parent_indices=NULL,
                                       const std::vector<uint>* sibling_indices=NULL);

    void fill_pc_fieldsplit_(const std::string &optionpath, PC &pc,  // fill the information about a fieldsplit pc
                             const std::string prefix, 
                             const uint &parent_offset, 
                             const std::vector<uint>* parent_indices);

    void fill_nullspace_(const std::string &optionpath, 
                         MatNullSpace &SP,
                         const uint &parent_offset, 
                         const std::vector<uint>* parent_indices);   // fill a petsc null space object

    void fill_constraints_();                                        // fill constraints on snes vi

    void fill_bound_(const std::string &optionpath, 
                         PETScVector_ptr &bound, const double &background_value);   // fill petsc vectors containing the bounds for snes vi

    void field_restrictions_(const std::string &optionpath,
                             std::vector<int>* &components,
                             std::vector<int>* &region_ids,
                             std::vector<int>* &boundary_ids,
                             dolfin::Expression* &value_exp,
                             double* &value_const,
                             const std::vector<std::size_t> &fieldshape,
                             const bool &fieldsymmetric);            // set up the restrictions on an IS by field

    void destroy_field_restrictions_(std::vector<int>* &components,  // destroy the objects describing any restrictions on an IS by field
                                     std::vector<int>* &region_ids,
                                     std::vector<int>* &boundary_ids,
                                     dolfin::Expression* &value_exp,
                                     double* &value_const);

    void initialize_tensors_();                                      // fill the tensor data structures of the solver bucket

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                           // empty the derived and base class data structures

    
  };
 
  typedef std::shared_ptr< SpudSolverBucket > SpudSolverBucket_ptr;// define a (boost shared) pointer to a spud solver bucket class

}
#endif
