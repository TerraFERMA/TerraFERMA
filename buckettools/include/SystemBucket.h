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
#include "ReferencePoints.h"
#include "PythonPeriodicMap.h"
#include <dolfin.h>

namespace buckettools
{

  class Bucket;                                                      // predeclaration
  typedef boost::shared_ptr< Bucket > Bucket_ptr;                    // so we can predeclare a pointer to it
  class SystemBucket;                                                // predeclare the class itself
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  
  enum functionbucket_type { SOLVE_START, SOLVE_TIMELOOP, SOLVE_DIAGNOSTICS, SOLVE_NEVER };

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

    void solve();                                                    // solve the solvers in this system (in order)

    void update();                                                   // update the functions in this system at the end of a timestep

    void update_nonlinear();                                         // update the potentially nonlinear functions in this system

    void update_timedependent();                                     // update the potentially time dependent functions in this system

    const double maxchange();                                        // return the maximum change in the (requested) fields over the last timestep

    void updatechange();                                             // update the change function

    void resetcalculated();                                          // reset the calculated booleans

    void cap_values();                                               // cap the values of the fields in this system

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    virtual void copy_diagnostics(SystemBucket_ptr &system, 
                            Bucket_ptr &bucket) const;               // copy the data necessary for the diagnostics data file(s)

    void initialize_diagnostics() const;                             // initialize any diagnostic output from this system

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

    const PythonPeriodicMap_ptr periodicmap() const                  // return a (boost shared) pointer to the system periodic map
    { return periodicmap_; }

    const std::vector<std::size_t> masterids() const
    { return masterids_; }

    const std::vector<std::size_t> slaveids() const
    { return slaveids_; }

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

    const Function_ptr changefunction() const                        // return a (boost shared) pointer to the change in the system
    { return changefunction_; }                                      // function over a timestep

    const Function_ptr residualfunction() const                      // return a (boost shared) pointer to the residual in the system
    { return residualfunction_; }                                    // function

    const Function_ptr snesupdatefunction() const                    // return a (boost shared) pointer to the snes update of the system
    { return snesupdatefunction_; }                                  // function

    const Expression_ptr icexpression() const                        // return a constant (boost shared) pointer to the initial
    { return icexpression_; }                                        // condition expression for this system

    const File_ptr icfile() const                                    // return a constant (boost shared) pointer to the initial
    { return icfile_; }                                              // condition file for this system

    File_ptr& icfile()                                                // return a (boost shared) pointer to the initial
    { return icfile_; }                                              // condition file for this system

    Bucket* bucket()                                                 // return a pointer to the parent bucket
    { return bucket_; }

    const Bucket* bucket() const                                     // return a constant pointer to the parent bucket
    { return bucket_; }

    const int solve_location() const                                 // return an integer describing where this system is solved
    { return solve_location_; }

    const bool solved() const                                        // return a boolean indicating if this system has been solved
    { return *solved_; }                                              // for or not

    const PETScVector_ptr residual_vector() const;                   // return the residual of the last solver in the system

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

    int_FunctionBucket_it orderedfields_begin();                     // return an iterator to the beginning of the ordered fields

    int_FunctionBucket_const_it orderedfields_begin() const;         // return a constant iterator to the beginning of the ordered fields

    int_FunctionBucket_it orderedfields_end();                       // return an iterator to the end of the ordered fields

    int_FunctionBucket_const_it orderedfields_end() const;           // return a constant iterator to the end of the ordered fields

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

    int_FunctionBucket_it orderedcoeffs_begin();                     // return an iterator to the beginning of the ordered coeffs

    int_FunctionBucket_const_it orderedcoeffs_begin() const;         // return a constant iterator to the beginning of the ordered coeffs

    int_FunctionBucket_it orderedcoeffs_end();                       // return an iterator to the end of the ordered coeffs

    int_FunctionBucket_const_it orderedcoeffs_end() const;           // return a constant iterator to the end of the ordered coeffs

    //***************************************************************|***********************************************************//
    // Solver bucket data access
    //***************************************************************|***********************************************************//

    void register_solver(SolverBucket_ptr solver, 
                                           const std::string &name); // register a solver bucket with the given name

    SolverBucket_ptr fetch_solver(const std::string &name);          // return a (boost shared) pointer to a solver with the
                                                                     // given name

    const SolverBucket_ptr fetch_solver(const std::string &name)     // return a constant (boost shared) pointer to a solver with the
                                                          const;     // given name

    SolverBucket_it solvers_begin();                                 // return an iterator to the beginning of the solver buckets

    SolverBucket_const_it solvers_begin() const;                     // return a constant iterator to the beginning of the solver
                                                                     // buckets

    SolverBucket_it solvers_end();                                   // return an iterator to the end of the solver buckets

    SolverBucket_const_it solvers_end() const;                       // return a constant iterator to the end of the solver buckets

    int_SolverBucket_it orderedsolvers_begin();                      // return an iterator to the beginning of the ordered solver
                                                                     // buckets

    int_SolverBucket_const_it orderedsolvers_begin() const;          // return a constant iterator to the beginning of the ordered
                                                                     // solver buckets

    int_SolverBucket_it orderedsolvers_end();                        // return an iterator to the end of the ordered solver buckets

    int_SolverBucket_const_it orderedsolvers_end() const;            // return a constant iterator to the end of the ordered solver
                                                                     // buckets

    //***************************************************************|***********************************************************//
    // BC data access
    //***************************************************************|***********************************************************//

    std::vector< const dolfin::DirichletBC* >::iterator dirichletbcs_begin();// return an iterator to the beginning of the system bcs

    std::vector< const dolfin::DirichletBC* >::const_iterator dirichletbcs_begin()// return a constant iterator to the beginning of the system
                                                          const;     // bcs

    std::vector< const dolfin::DirichletBC* >::iterator dirichletbcs_end();  // return an iterator to the end of the system bcs

    std::vector< const dolfin::DirichletBC* >::const_iterator dirichletbcs_end()// return a constant iterator to the end of the system bcs
                                                          const;

    const std::vector< const dolfin::DirichletBC* > dirichletbcs() const     // return a constant vector of system bcs
    { return dirichletbcs_; }
    
    //***************************************************************|***********************************************************//
    // Reference point data access
    //***************************************************************|***********************************************************//

    std::vector<ReferencePoints_ptr>::iterator points_begin();        // return an iterator to the beginning of the system reference
                                                                     // points

    std::vector<ReferencePoints_ptr>::const_iterator points_begin()   // return a constant iterator to the beginning of the system
                                                          const;     // reference points

    std::vector<ReferencePoints_ptr>::iterator points_end();          // return an iterator to the end of the system reference points

    std::vector<ReferencePoints_ptr>::const_iterator points_end()     // return a constant iterator to the end of the system
                                                          const;     // reference points

    const std::vector< ReferencePoints_ptr > points() const           // return a constant vector of system reference points
    { return points_; }
    
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
    
    virtual const std::string str() const                            // return a string describing the contents of the system
    { return str(0); }

    virtual const std::string str(int indent) const;                 // return an indented string describing the contents of the
                                                                     // system

    virtual const std::string fields_str() const                     // return a string describing the fields in the system
    { return fields_str(0); }

    virtual const std::string fields_str(const int &indent) const;   // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string coeffs_str() const                     // return a string describing the coefficients in the system
    { return coeffs_str(0); }

    virtual const std::string coeffs_str(const int &indent) const;   // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string solvers_str() const                    // return a string describing the solver buckets in the system
    { return coeffs_str(0); }

    virtual const std::string solvers_str(const int &indent) const;  // return an indented string describing the solver buckets in
                                                                     // the system

    void checkpoint();                                               // checkpoint the system

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_all_coeffs_();                                       // attach all fields and coefficients to forms and functionals

    void attach_function_coeffs_(FunctionBucket_it f_begin,          // attach specific fields or coefficients to functionals
                                          FunctionBucket_it f_end);

    void attach_solver_coeffs_(SolverBucket_it s_begin,              // attach specific fields or coefficients to solver forms
                                          SolverBucket_it s_end);

    void collect_ics_(const uint &component, const std::size_t &maxrank, // collect the field initial conditions into an initial 
                      const std::map< std::size_t, Expression_ptr >  // condition expression
                                                  &icexpressions);

    void apply_ic_();                                                // apply the initial conditions to the system function

    void apply_dirichletbc_();                                       // apply the Dirichlet bcs to the system function

    void apply_referencepoints_();                                   // apply the reference points to the system function

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the system name

    std::string uflsymbol_;                                          // the system ufl symbol

    Bucket* bucket_;                                                 // a pointer to the parent bucket

    Mesh_ptr mesh_;                                                  // a (boost shared) pointer to the system mesh

    MeshFunction_size_t_ptr celldomains_, facetdomains_;             // a (boost shared) pointer to the system mesh functions

    FunctionSpace_ptr functionspace_;                                // a (boost shared) pointer to the system functionspace

    Function_ptr function_, oldfunction_, iteratedfunction_;         // (boost shared) pointers to the system functions at different
                                                                     // time levels (old, iterated - most up to date -, base)

    Expression_ptr icexpression_;                                    // (boost shared) pointer to an expression describing the initial condition

    File_ptr icfile_;                                                // (boost shared) pointer to a file containing a checkpointed ic

    int solve_location_;                                             // when this system will be solved

    Function_ptr changefunction_;                                    // (boost shared) pointer to the change between timesteps

    bool_ptr change_calculated_;                                     // indicate if the change has been recalculated recently

    bool_ptr solved_;                                                // indicate if the system has been solved this timestep

    Function_ptr residualfunction_;                                  // (boost shared) pointer to the residual of the system

    Function_ptr snesupdatefunction_;                                // (boost shared) pointer to the snes update of the system

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, FunctionBucket_ptr > fields_;             // a map from field names to (boost shared) pointers to fields
    
    std::map< int, FunctionBucket_ptr > orderedfields_;              // an ordered (user defined)  map to (boost shared) pointers to fields
    
    std::map< std::string, FunctionBucket_ptr > coeffs_;             // a map from coefficient names to (boost shared) pointers to
                                                                     // coefficients

    std::map< int, FunctionBucket_ptr > orderedcoeffs_;              // an ordered (user defined)  map to (boost shared) pointers to coeffs
    
    std::map< std::string, SolverBucket_ptr > solvers_;              // a map from solver bucket names to (boost shared) pointers to
                                                                     // solver buckets

    std::map< int, SolverBucket_ptr > orderedsolvers_;               // an ordered (user defined) map to
                                                                     // (boost shared) pointers to solver buckets

    std::vector< const dolfin::DirichletBC* > dirichletbcs_;         // a vector of (boost shared) poitners to the dirichlet bcs

    PythonPeriodicMap_ptr periodicmap_;                              // periodic map for (a single) periodic bc

    std::vector<std::size_t> masterids_, slaveids_;                  // master and slave ids for periodic bc

    std::vector< ReferencePoints_ptr > points_;                      // a vector of (boost shared) poitners to reference points

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    virtual void checkpoint_options_();                              // checkpoint the options system for the systembucket

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the data structures in this system

  };

  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // define a (boost shared) pointer to the system class type

}
#endif
