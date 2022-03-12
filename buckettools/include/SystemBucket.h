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


#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include "FunctionalBucket.h"
#include <dolfin.h>

namespace buckettools
{

  class Bucket;                                                      // predeclaration
  typedef std::shared_ptr< Bucket > Bucket_ptr;                    // so we can predeclare a pointer to it
  class SystemBucket;                                                // predeclare the class itself
  typedef std::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  
  //*****************************************************************|************************************************************//
  // SystemBucket class:
  //
  // The SystemBucket class describes a functionspace and a set of solvers that act on the fields
  // contained in that (potentially mixed) functionspace.
  // This base class describes the basic data structures while derived classes may be defined
  // that allow it to be linked to an options system.
  //*****************************************************************|************************************************************//
  class SystemBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SystemBucket();                                                  // default constructor

    SystemBucket(Bucket* bucket);                                    // specific constructor
    
    virtual ~SystemBucket();                                         // default destructor

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    void evaluate_initial_fields();                                  // evaluate the initial values of the fields

    bool solve(const int &location, const bool force=true);          // solve the solvers in this system (in order)

    bool solve(const std::vector<int> &locations=std::vector<int>(), 
                                 const bool force=true);             // solve the solvers in this system matching any locations

    void update();                                                   // update the functions in this system at the end of a timestep

    void update_nonlinear();                                         // update the potentially nonlinear functions in this system

    void update_iterated(const double relax);                        // update the iterated vectors in this system

    void update_timedependent();                                     // update the potentially time dependent functions in this system

    const double maxchange();                                        // return the maximum change in the (requested) fields over the last timestep

    void updatechange();                                             // update the change function

    void resetcalculated();                                          // reset the calculated booleans

    void postprocess_values();                                       // cap the values of the fields in this system

    double residual_norm(const int &location);                       // return the norm of the residual of the last solver that meets the location requirement

    double residual_norm(const std::vector<int> &locations=
                                                 std::vector<int>());// return the norm of the residual of the last solver that meets the location requirement

    const std::vector<int> solve_locations() const;                  // return a std::vector indicating the solve locations of the solvers

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void initialize_forms();                                         // attach all fields and coefficients to forms and functionals

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const Function_ptr function_ptr(const double_ptr time) const;    // return a (boost shared) pointer to the system function

    const std::string name() const                                   // return the name of this system
    { return name_; }

    const std::string uflsymbol() const                              // return the system ufl symbol
    { return uflsymbol_; }

    const Mesh_ptr mesh() const                                      // return a (boost shared) pointer to the system mesh
    { return mesh_; }

    const MeshFunction_size_t_ptr celldomains() const                // return a (boost shared) pointer to the system mesh function
    { return celldomains_; }

    const MeshFunction_size_t_ptr facetdomains() const               // return a (boost shared) pointer to the system mesh function
    { return facetdomains_; }

    const FunctionSpace_ptr functionspace() const                    // return a (boost shared) pointer to the system functionspace
    { return functionspace_; }

    const Function_ptr function() const                              // return a (boost shared) pointer to the system function
    { return function_; }

    const Function_ptr oldfunction() const                           // return a (boost shared) pointer to the old system function
    { return oldfunction_; }

    const Function_ptr iteratedfunction() const                      // return a (boost shared) pointer to the iterated system
    { return iteratedfunction_; }                                    // function

    const Function_ptr olditeratedfunction() const                   // return a (boost shared) pointer to the "old" iterated system
    { return olditeratedfunction_; }                                 // function from the previous iteration

    const Function_ptr changefunction() const                        // return a (boost shared) pointer to the change in the system
    { return changefunction_; }                                      // function over a timestep

    const Function_ptr residualfunction() const                      // return a (boost shared) pointer to the residual in the system
    { return residualfunction_; }                                    // function

    const Function_ptr snesupdatefunction() const                    // return a (boost shared) pointer to the snes update of the system
    { return snesupdatefunction_; }                                  // function

    Bucket* bucket()                                                 // return a pointer to the parent bucket
    { return bucket_; }

    const Bucket* bucket() const                                     // return a constant pointer to the parent bucket
    { return bucket_; }

    //***************************************************************|***********************************************************//
    // Field data access
    //***************************************************************|***********************************************************//

    void register_field(FunctionBucket_ptr field, 
                                          const std::string &name); // register a field (subfunction) with the given name

    FunctionBucket_ptr fetch_field(const std::string &name);         // return a (boost shared) pointer to a field with the given
                                                                     // name

    const FunctionBucket_ptr fetch_field(const std::string &name)    // return a constant (boost shared) pointer to a field with the
                                                          const;     // given name

    FunctionBucket_it fields_begin();                                // return an iterator to the beginning of the fields

    FunctionBucket_const_it fields_begin() const;                    // return a constant iterator to the beginning of the fields

    FunctionBucket_it fields_end();                                  // return an iterator to the end of the fields

    FunctionBucket_const_it fields_end() const;                      // return a constant iterator to the end of the fields

    const int fields_size() const;                                   // return the number of fields

    //***************************************************************|***********************************************************//
    // Coefficient data access
    //***************************************************************|***********************************************************//

    void register_coeff(FunctionBucket_ptr coeff, 
                                          const std::string &name); // register a coefficient with the given name

    FunctionBucket_ptr fetch_coeff(const std::string &name);         // return a (boost shared) pointer to a coefficient with the
                                                                     // given name

    const FunctionBucket_ptr fetch_coeff(const std::string &name)    // return a constant (boost shared) pointer to a coefficient with the
                                                          const;     // given name

    FunctionBucket_it coeffs_begin();                                // return an iterator to the beginning of the coefficients

    FunctionBucket_const_it coeffs_begin() const;                    // return a constant iterator to the beginning of the
                                                                     // coefficients

    FunctionBucket_it coeffs_end();                                  // return an iterator to the end of the coefficients

    FunctionBucket_const_it coeffs_end() const;                      // return a constant iterator to the end of the coefficients

    //***************************************************************|***********************************************************//
    // Solver bucket data access
    //***************************************************************|***********************************************************//

    void register_solver(SolverBucket_ptr solver, 
                                          const std::string &name);        // register a solver bucket with the given name

    SolverBucket_ptr fetch_solver(const std::string &name);          // return a (boost shared) pointer to a solver with the
                                                                     // given name

    const SolverBucket_ptr fetch_solver(const std::string &name)     // return a constant (boost shared) pointer to a solver with the
                                                          const;     // given name

    SolverBucket_it solvers_begin();                                 // return an iterator to the beginning of the solver buckets

    SolverBucket_const_it solvers_begin() const;                     // return a constant iterator to the beginning of the solver
                                                                     // buckets

    SolverBucket_it solvers_end();                                   // return an iterator to the end of the solver buckets

    SolverBucket_const_it solvers_end() const;                       // return a constant iterator to the end of the solver buckets

    //***************************************************************|***********************************************************//
    // Functional data access
    //***************************************************************|***********************************************************//

    void register_functional(FunctionalBucket_ptr functional, 
                                           const std::string &name); // register a functional with the given name in this function

    FunctionalBucket_ptr fetch_functional(const std::string &name);  // return a (boost shared) pointer to a functional with the given name

    const FunctionalBucket_ptr fetch_functional(const std::string &name) 
                                                             const;  // return a constant (boost shared) pointer to a functional with the given name

    FunctionalBucket_it functionals_begin();                         // return an iterator to the beginning of the functionals of this function

    FunctionalBucket_const_it functionals_begin() const;             // return a constant iterator to the beginning of the functionals of this function

    FunctionalBucket_it functionals_end();                           // return an iterator to the end of the functionals of this function

    FunctionalBucket_const_it functionals_end() const;               // return a constant iterator to the end of the functionals of this function

    //***************************************************************|***********************************************************//
    // BC data access
    //***************************************************************|***********************************************************//

    std::vector< std::shared_ptr<const dolfin::DirichletBC> >::iterator bcs_begin();// return an iterator to the beginning of the system bcs

    std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator bcs_begin()// return a constant iterator to the beginning of the system
                                                          const;     // bcs

    std::vector< std::shared_ptr<const dolfin::DirichletBC> >::iterator bcs_end();  // return an iterator to the end of the system bcs

    std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator bcs_end()// return a constant iterator to the end of the system bcs
                                                          const;

    const std::vector< std::shared_ptr<const dolfin::DirichletBC> > bcs() const     // return a constant vector of system bcs
    { return bcs_; }
    
    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const bool include_in_visualization() const;                     // return a boolean indicating if this system has fields to 
                                                                     // be included in diagnostic output
    
    const bool include_in_statistics() const;                        // return a boolean indicating if this system has fields to 
                                                                     // be included in diagnostic output
    
    const bool include_in_steadystate() const;                       // return a boolean indicating if this system has fields to 
                                                                     // be included in steadystate output
    
    const bool include_in_detectors() const;                         // return a boolean indicating if this system has fields to 
                                                                     // be included in steadystate output
    
    virtual const std::string str(int indent=0) const;               // return an indented string describing the contents of the
                                                                     // system

    virtual const std::string fields_str(const int &indent=0) const; // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string coeffs_str(const int &indent=0) const; // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string solvers_str(const int &indent=0) const;// return an indented string describing the solver buckets in
                                                                     // the system

    virtual const std::string functionals_str(const int &indent=0) const;// return an indented string describing the functionals 
                                                                     // of the system

    void output();                                                   // output vis from the system

    void checkpoint(const double_ptr time);                          // checkpoint the system

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void apply_ic_();                                                // apply the initial conditions to the system function

    void apply_bcs_();                                               // apply the Dirichlet bcs to the system function

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the system name

    std::string uflsymbol_;                                          // the system ufl symbol

    Bucket* bucket_;                                                 // a pointer to the parent bucket

    Mesh_ptr mesh_;                                                  // a (boost shared) pointer to the system mesh

    MeshFunction_size_t_ptr celldomains_, facetdomains_;             // a (boost shared) pointer to the system mesh functions

    FunctionSpace_ptr functionspace_;                                // a (boost shared) pointer to the system functionspace

    Function_ptr function_, oldfunction_, iteratedfunction_,         // (boost shared) pointers to the system functions at different
                 olditeratedfunction_;                               // time levels (old, iterated - most up to date -, base,
                                                                     // olditerated - previous iteration)

    Function_ptr changefunction_;                                    // (boost shared) pointer to the change between timesteps

    bool_ptr change_calculated_;                                     // indicate if the change has been recalculated recently

    Function_ptr residualfunction_;                                  // (boost shared) pointer to the residual of the system

    Function_ptr snesupdatefunction_;                                // (boost shared) pointer to the snes update of the system

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    ordered_map<const std::string, FunctionBucket_ptr> fields_;      // a map from field names to (boost shared) pointers to fields
    
    ordered_map<const std::string, FunctionBucket_ptr> coeffs_;      // a map from coefficient names to (boost shared) pointers to
                                                                     // coefficients

    ordered_map<const std::string, SolverBucket_ptr> solvers_;       // a map from solver bucket names to (boost shared) pointers to
                                                                     // solver buckets

    ordered_map<const std::string, FunctionalBucket_ptr> functionals_;// map from functional names to form (boost shared) pointers

    std::vector< std::shared_ptr<const dolfin::DirichletBC> > bcs_;                  // a vector of (boost shared) poitners to the dirichlet bcs

  };

  typedef std::shared_ptr< SystemBucket > SystemBucket_ptr;        // define a (boost shared) pointer to the system class type

}
#endif
