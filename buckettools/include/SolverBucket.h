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


#ifndef __SOLVERBUCKET_H
#define __SOLVERBUCKET_H

#include "BoostTypes.h"
#include "BucketPETScBase.h"
#include "ConvergenceFile.h"
#include "KSPConvergenceFile.h"
#include "GenericSolverBucket.h"
#include "SystemBucket.h"
#include <dolfin.h>
#include "petscsnes.h"

namespace buckettools
{

  class SystemBucket;                                                // predeclare
  typedef std::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  class FunctionBucket;                                              // predeclare the class itself
  typedef std::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // so we can predeclare a pointer to it

  //*****************************************************************|************************************************************//
  // SolverBucket class:
  //
  // The SolverBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  //*****************************************************************|************************************************************//
  class SolverBucket : public GenericSolverBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SolverBucket();                                                  // default constructor

    SolverBucket(SystemBucket* system);                              // specific constructor (with parent system)
    
    virtual ~SolverBucket();                                         // default destructor

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    bool solve();                                                    // run the nonlinear solver described by this class

    double residual_norm();                                          // return the norm of the residual (which will be reassembled)

    void resetcalculated();                                          // update this solver at the end of a timestep

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_form_coeffs();                                       // attach coefficients to the forms in this solver

    void create_nullspace();                                         // take any stored nullspace vectors and convert them into a
                                                                     // PETSc null space object

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string type() const                                   // return a string describing the solver type
    { return type_; }

    const SNES snes() const                                          // return the snes object being used
    { return snes_; }

    const Form_ptr bilinear_form() const                             // return a (boost shared) pointer to the bilinear form
    { return bilinear_; }

    const bool ident_zeros() const                                   // return true if the bilinear form needs to be idented
    { return ident_zeros_; }

    const Form_ptr bilinearpc_form() const                           // return a (boost shared) pointer to the bilinear pc form
    { return bilinearpc_; }

    const bool ident_zeros_pc() const                                // return true if the bilinear pc form needs to be idented
    { return ident_zeros_pc_; }

    const Form_ptr linear_form() const                               // return a (boost shared) pointer to the linear form
    { return linear_; }

    const PETScVector_ptr residual_vector() const                    // return the residual of this solver
    { return res_; }

    SystemBucket* system()                                           // return a pointer to the parent system
    { return system_; }

    const SystemBucket* system() const                               // return a pointer to the parent system
    { return system_; }

    const int iteration_count() const;                               // return the number of iterations taken (may not be accurate!)

    void iteration_count(const int &it);                             // set the number of iterations taken

    const bool visualization_monitor() const;                        // return true if we're using a visualization monitor
    
    const bool kspvisualization_monitor() const;                     // return true if we're using a visualization monitor

    std::string visualization_basename();                            // return a basename for any visualization monitors
    
    const bool monitor_norms() const                                 // return true if norms should be monitored in nonlinear iterations
    { return monitornorms_; }

    MatNullSpace nullspace()
    { return sp_; }

    const int solve_location() const                                 // return an integer describing where this solver is applied
    { return solve_location_; }

    std::string current_systemssolver() const                        // return a flag indicating which systems solver we're
    { return current_systemssolver_; }                               // currently using (or assuming we're using)

    void set_current_systemssolver(const std::string systemssolvername)
    { current_systemssolver_ = systemssolvername; }

    void reset_current_systemssolver()
    { current_systemssolver_ = ""; }

    //***************************************************************|***********************************************************//
    // Form data access
    //***************************************************************|***********************************************************//

    void register_form(Form_ptr form, const std::string &name);      // register a form in the solver

    bool contains_form(const std::string &name);                     // return a boolean, true if the names form exists in the
                                                                     // solver, false otherwise

    Form_ptr fetch_form(const std::string &name);                    // fetch the named form

    Form_it forms_begin();                                           // return an iterator to the beginning of the forms

    Form_const_it forms_begin() const;                               // return a constant iterator to the beginning of the forms

    Form_it forms_end();                                             // return an iterator to the end of the forms

    Form_const_it forms_end() const;                                 // return a constant iterator to the end of the forms

    //***************************************************************|***********************************************************//
    // Solver form data access
    //***************************************************************|***********************************************************//

    void register_solverform(Form_ptr form, const std::string &name);// register a form in the solver

    Form_it solverforms_begin();                                     // return an iterator to the beginning of the forms

    Form_const_it solverforms_begin() const;                         // return a constant iterator to the beginning of the forms

    Form_it solverforms_end();                                       // return an iterator to the end of the forms

    Form_const_it solverforms_end() const;                           // return a constant iterator to the end of the forms

    PETScMatrix_ptr fetch_solvermatrix(const std::string &name);     // fetch the named solver matrix

    bool solverident_zeros(const std::string &name);                 // ident zero the named solver matrix

    IS fetch_solverindexset(const std::string &name);                // fetch the named ident zeros

    Mat fetch_solversubmatrix(const std::string &name);              // fetch the named solver submatrix

    //***************************************************************************************************************************//
    // Parent nonlinear systems solver data access
    //***************************************************************************************************************************//

    void register_systemssolver(SystemsSolverBucket* solver, const std::string &name);//register a systems solver
    
    SystemsSolverBucket* fetch_systemssolver(const std::string &name);// fetch the named systemssolver

    p_SystemsSolverBucket_it systemssolvers_begin();                 // return an iterator to the beginning of the systemssolvers

    p_SystemsSolverBucket_const_it systemssolvers_begin() const;     // return a constant iterator to the beginning of the systemssolvers

    p_SystemsSolverBucket_it systemssolvers_end();                   // return an iterator to the end of the systemssolvers

    p_SystemsSolverBucket_const_it systemssolvers_end() const;       // return a constant iterator to the end of the systemssolvers

    void register_convergencefile(ConvergenceFile_ptr convfile, const std::string &name="");//register a convergence file
    
    ConvergenceFile_ptr convergencefile();                           // get the default convergence file

    ConvergenceFile_ptr fetch_convergencefile(const std::string &name);// fetch the named convergence file

    ConvergenceFile_it convergencefiles_begin();                     // return an iterator to the beginning of the convergencefiles

    ConvergenceFile_const_it convergencefiles_begin() const;         // return a constant iterator to the beginning of the convergencefiles

    ConvergenceFile_it convergencefiles_end();                       // return an iterator to the end of the convergencefiles

    ConvergenceFile_const_it convergencefiles_end() const;           // return a constant iterator to the end of the convergencefiles

    void register_kspconvergencefile(KSPConvergenceFile_ptr kspconvfile, const std::string &name="");//register a ksp convergence file

    KSPConvergenceFile_ptr kspconvergencefile();                     // get the default ksp convergence file
    
    KSPConvergenceFile_ptr fetch_kspconvergencefile(const std::string &name);// fetch the named ksp convergence file

    KSPConvergenceFile_it kspconvergencefiles_begin();               // return an iterator to the beginning of the ksp convergencefiles

    KSPConvergenceFile_const_it kspconvergencefiles_begin() const;   // return a constant iterator to the beginning of the ksp convergencefiles

    KSPConvergenceFile_it kspconvergencefiles_end();                 // return an iterator to the end of the ksp convergencefiles

    KSPConvergenceFile_const_it kspconvergencefiles_end() const;     // return a constant iterator to the end of the ksp convergencefiles

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    virtual const std::string str(int indent=0) const;               // return an indented string describing the contents of this
                                                                     // solver

    virtual const std::string forms_str(const int &indent=0) const;  // return an indented string describing the forms in this
                                                                     // solver

    void checkpoint();                                               // checkpoint the solverbucket

    const bool include_in_visualization() const                      // return a boolean indicating if this system has fields to 
    { return (*system_).include_in_visualization(); }                // be included in diagnostic output
    
    const bool include_in_statistics() const                         // return a boolean indicating if this system has fields to 
    { return (*system_).include_in_statistics(); }                   // be included in diagnostic output
    
    const bool include_in_steadystate() const                        // return a boolean indicating if this system has fields to 
    { return (*system_).include_in_steadystate(); }                  // be included in steadystate output
    
    const bool include_in_detectors() const                          // return a boolean indicating if this system has fields to 
    { return (*system_).include_in_detectors(); }                    // be included in steadystate output
    

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    SNES snes_;                                                      // a petsc snes object

    SNESCtx ctx_;                                                    // a snes context type (defined in BucketPETScBase.h)

    KSP ksp_;                                                        // a petsc ksp type

    Form_ptr linear_, bilinear_, bilinearpc_, residual_;             // (boost shared) pointers to forms

    ordered_map<const std::string, Form_ptr> solverforms_;                  // (boost shared) pointers to forms under linear solver 

    PETScMatrix_ptr matrix_, matrixpc_;                              // dolfin petsc matrix types

    std::map< std::string, PETScMatrix_ptr > solvermatrices_;        // dolfin petsc matrices for solver matrices

    std::map< std::string, IS > solverindexsets_;                    // (boost shared) pointers to the indexsets defining the solver submatrices

    std::map< std::string, Mat > solversubmatrices_;                 // (boost shared) pointers to the sub petsc matrices

    PETScVector_ptr rhs_, res_, work_;                               // dolfin petsc vector types

    double rtol_, atol_, stol_;                                      // nonlinear solver tolerances

    int minits_, maxits_, maxfes_;                                   // nonlinear solver iteration counts

    int_ptr iteration_count_;                                        // nonlinear iterations taken (may not be accurate!)

    bool ident_zeros_, ident_zeros_pc_;                              // replace zero rows with the identity (matrix and pc)

    std::map< std::string, bool > solverident_zeros_;                // replace zero rows with the identity (solver matrices)

    bool ignore_failures_;                                           // ignore solver failures

    std::string type_;                                               // solver type (string)

    SystemBucket* system_;                                           // parent system

    CustomMonitorCtx snesmctx_, kspmctx_;                            // monitor contexts

    bool_ptr visualizationmonitor_, kspvisualizationmonitor_;        // visualization monitors

    bool monitornorms_;                                              // monitor the norms in nonlinear iterations

    int solve_location_;                                             // when this solver will be applied

    std::string current_systemssolver_;                              // an indicator of which systemssolver we're currently using

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    ordered_map<const std::string, Form_ptr> forms_;                 // a map from the form names to the form pointers

    ordered_map<const std::string, SystemsSolverBucket*> systemssolvers_;// a map of the systems solvers using this solver

    ordered_map<const std::string, ConvergenceFile_ptr> convergencefiles_;  // a map from systems solvers names to convergence files

    ordered_map<const std::string, KSPConvergenceFile_ptr> kspconvergencefiles_;// a map from systems solvers names to ksp convergence files

    std::vector< PETScVector_ptr > nullspacevectors_;                // a vector of null space vectors to be removed from the rhs
                                                                     // after assembly

    MatNullSpace sp_;                                                // PETSc matnullspace object

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void orthonormalize_petsc_vecs_(Vec vecs[], PetscInt n);         // orthonormalize an array of PETSc Vecs

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    virtual void checkpoint_options_();                              // checkpoint the options system for the solverbucket

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    //***************************************************************|***********************************************************//
    // Solver convergence checking
    //***************************************************************|***********************************************************//

    void snes_check_convergence_();                                  // check snes convergence

    void ksp_check_convergence_(KSP &ksp, int indent);               // check ksp convergence

    void ksp_check_convergence_(KSP &ksp)                            // check ksp convergence (no indent)
    { ksp_check_convergence_(ksp, 0); }

  };

  typedef std::shared_ptr< SolverBucket > SolverBucket_ptr;        // define a (boost shared) pointer to the solver class type

}
#endif
