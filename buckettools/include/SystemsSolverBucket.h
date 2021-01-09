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


#ifndef __SYSTEMSSOLVERBUCKET_H
#define __SYSTEMSSOLVERBUCKET_H

#include "GenericSolverBucket.h"
#include "BoostTypes.h"
#include "SystemsConvergenceFile.h"

namespace buckettools
{
  
  //class SolverBucket;                                                // predeclare
  //typedef std::shared_ptr< SolverBucket > SolverBucket_ptr;

  //*****************************************************************|************************************************************//
  // SystemsSolverBucket class:
  //
  //*****************************************************************|************************************************************//
  class SystemsSolverBucket : public GenericSolverBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SystemsSolverBucket();                                           // default constructor

    SystemsSolverBucket(const int &solve_location, Bucket* bucket, 
                        SystemsSolverBucket* systemssolver=NULL);    // specific constructor (with parent bucket and solver)
    
    virtual ~SystemsSolverBucket();                                  // default destructor


    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill the SystemsSolverBucket

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    bool solve();                                                    // run the nonlinear solver described by this class

    bool solve_diagnostics(const bool &write_vis,                    // a special interface to solve for diagnostics
                           const bool &write_stat, 
                           const bool &write_steady, 
                           const bool &write_det);

    double residual_norm();                                          // return the norm of the residual (which will be reassembled)

    const std::vector<SolverBucket_ptr>& residualsolvers() const
    { return residualsolvers_; }

    void resetcalculated();                                          // update this solver at the end of a timestep

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string type() const                                   // return a string containing the solver type
    { return type_; }

    Bucket* bucket()                                                 // return a pointer to the parent bucket
    { return bucket_; }

    const Bucket* bucket() const                                     // return a pointer to the parent bucket
    { return bucket_; }

    SystemsSolverBucket* systemssolver()                             // return a pointer to the parent systems solver
    { return systemssolver_; }

    const SystemsSolverBucket* systemssolver() const                 // return a pointer to the parent systems solver
    { return systemssolver_; }

    const int iteration_count() const;                               // return the number of iterations taken (may not be accurate!)

    void iteration_count(const int &it);                             // set the number of iterations taken

    const int solve_location() const                                 // return an integer describing where this solver is applied
    { return solve_location_; }

    const bool iterative() const                                     // return an indicator of whether we are iterating (based on
    { return rtol_; }                                                // rtol_'s allocation)

    const bool coupled() const                                       // return an indicator of whether the systems are coupled
    { return coupled_; }

    std::string unique_name() const;

    std::string iterations_str() const;                              //utility function to get the iterations as string

    //***************************************************************|***********************************************************//
    // Sub solver data access
    //***************************************************************|***********************************************************//

    void register_solver(GenericSolverBucket_ptr solver, const std::string &name);// register a sub solver in the solver

    GenericSolverBucket_ptr fetch_solver(const std::string &name);   // fetch the named solver

    const GenericSolverBucket_ptr fetch_solver(const std::string &name) const;// fetch the named solver

    GenericSolverBucket_it solvers_begin();                          // return an iterator to the beginning of the solvers

    GenericSolverBucket_const_it solvers_begin() const;              // return a constant iterator to the beginning of the solvers

    GenericSolverBucket_it solvers_end();                            // return an iterator to the end of the solvers

    GenericSolverBucket_const_it solvers_end() const;                // return a constant iterator to the end of the solvers

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    virtual const std::string str(int indent=0) const;               // return an indented string describing the contents of this

    virtual const std::string solvers_str(const int &indent) const;  // return an indented string describing the sub solvers

    void checkpoint();                                               // checkpoint the systemssolverbucket

    const bool include_in_visualization() const;                     // return a boolean indicating if this systems solver has fields to 
                                                                     // be included in diagnostic output
    
    const bool include_in_statistics() const;                        // return a boolean indicating if this systems solver has fields to 
                                                                     // be included in diagnostic output
    
    const bool include_in_steadystate() const;                       // return a boolean indicating if this systems solver has fields to 
                                                                     // be included in steadystate output
    
    const bool include_in_detectors() const;                         // return a boolean indicating if this systems solver has fields to 
                                                                     // be included in steadystate output
    
  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill the base data of the systems solver bucket

    void fill_solvers_();                                            // fill the child solvers

    void fill_solverlists_();                                        // fill information about the solvers

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string type_;                                               // solver type (string)

    int solve_location_;                                             // when this solver will be applied

    Bucket* bucket_;                                                 // parent bucket

    SystemsSolverBucket* systemssolver_;                             // parent systems solver (if it exists)

    int_ptr iteration_count_;                                        // the number of iterations requested and the number of nonlinear 
                                                                     // iterations taken
    int minits_, maxits_;                                            // nonlinear system iteration counts

    double *rtol_, atol_;                                            // nonlinear system tolerances

    bool ignore_failures_;                                           // ignore solver failures

    bool coupled_;                                                   // are the solvers coupled?

    //***************************************************************|***********************************************************//
    // Diagnostics data
    //***************************************************************|***********************************************************//

    SystemsConvergenceFile_ptr convfile_;                            // nonlinear systems convergence file

    bool write_convvis_;                                             // write convergence visualization files every iteration

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    ordered_map<const std::string, GenericSolverBucket_ptr> solvers_;// a map from the solver names to the solver pointers

    std::vector<SolverBucket_ptr> residualsolvers_;                  // a list of the solvers used to caluclate the residual norm

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to members of this class

    //***************************************************************|***********************************************************//
    // Solver convergence checking
    //***************************************************************|***********************************************************//

    bool complete_iterating_(const double &aerror0);                 // indicate if nonlinear systems iterations are complete or not

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    virtual void checkpoint_options_();                              // checkpoint the options system for the systemssolverbucket

  };

  typedef std::shared_ptr< SystemsSolverBucket > SystemsSolverBucket_ptr;// define a boost shared ptr type for the class

}
#endif

