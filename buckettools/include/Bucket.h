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


#ifndef __BUCKET_H
#define __BUCKET_H

#include "GenericDetectors.h"
#include "BoostTypes.h"
#include "SystemBucket.h"
#include "StatisticsFile.h"
#include "SteadyStateFile.h"
#include "SolverBucket.h"
#include "DetectorsFile.h"
#include "SystemsConvergenceFile.h"
#include <dolfin.h>
#include <boost/timer/timer.hpp>

namespace buckettools
{
  
  class Bucket;                                                      // predeclare the class itself
  typedef std::shared_ptr< Bucket > Bucket_ptr;                    // so we can predeclare a pointer to it
  
  enum output_location { OUTPUT_START, OUTPUT_TIMELOOP, OUTPUT_END };

  enum checkpoint_location { CHECKPOINT_TIMELOOP, CHECKPOINT_END };

  //*****************************************************************|************************************************************//
  // Bucket class:
  //
  // A class that describes a collection of meshes, systems, ufl symbols and point diagnostics (detectors).
  // The base class contains data structures and mechanisms for accessing those structures.
  // Derived classes are designed for specific options systems and are capable of populating those data
  // structures and may be extended to incorporate the ability to run entire models.
  //*****************************************************************|************************************************************//
  class Bucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    Bucket();                                                        // default constructor

    Bucket(const std::string &name);                                 // specific constructor specifying the name
    
    virtual ~Bucket();                                               // default destructor

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    void run();                                                      // run the model described in the bucket
 
    void solve(const int &location);                                 // solve all the systems in the bucket (depending on location)

    void resetcalculated();                                          // update the calculated booleans in the systems in this bucket

    void update();                                                   // update the functions in the systems in this bucket

    void update_timestep();                                          // update the timestep

    void update_nonlinear();                                         // update the potentially nonlinear functions in the systems in this bucket

    void update_iterated();                                          // update the iterated vectors in the systems (allows relaxation in global nonlinear iterations)

    void update_timedependent();                                     // update the potentially timedependent functions in the systems in this bucket

    bool complete();                                                 // indicate if the simulation is complete or not

    bool complete_timestepping();                                    // indicate if timestepping is complete or not

    double residual_norm();                                          // return the l2 norm of the residual of the systems

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_coeffs(Form_it f_begin, Form_it f_end);              // attach coefficients to the selected forms

    void attach_coeffs(Form_ptr form);                               // attach coefficients to the selected form

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a string containing the name of the bucket
    { return name_; }
    
    const int dimension() const                                      // return the geometry dimension
    { return dimension_; }

    const int timestep_count() const;                                // return the timestep count

    const double start_time() const;                                 // return the current time

    const double_ptr start_time_ptr() const                          // return a (const std shared) pointer to the start time
    { return start_time_; }

    const double old_time() const;                                   // return the old time (from the previous timestep)

    const double_ptr old_time_ptr() const                            // return a (const std shared) pointer to the old time
    { return old_time_; }

    const double current_time() const;                               // return the current time

    const double_ptr current_time_ptr() const                        // return a (const std shared) pointer to the current time
    { return current_time_; }

    const double finish_time() const;                                // return the finish time

    const int number_timesteps() const;                              // return the number of timesteps after which the simulation will finish

    const double walltime_limit() const;                             // return the walltime limit

    const double timestep() const;                                   // return the timestep (as a double)

    const int iteration_count() const;                               // return the number of nonlinear iterations taken

    const double relaxation_parameter() const                        // return the relaxation parameter
    { return relax_; }

    const std::string output_basename() const                        // return the output base name
    { return output_basename_; }

    static const time_t* start_walltime()                            // return the start time
    { return &start_walltime_; }

    static const double elapsed_walltime()                           // return the elapsed wall time
    { return dolfin::MPI::max(MPI_COMM_WORLD, static_cast<double>(timer_.elapsed().wall)*1.e-9); }

    const int checkpoint_count() const;                              // return the checkpoint count

    //***************************************************************|***********************************************************//
    // Mesh data access
    //***************************************************************|***********************************************************//

    void register_mesh(Mesh_ptr mesh, const std::string &name,
                       MeshFunction_size_t_ptr celldomains  = NULL,
                       MeshFunction_size_t_ptr facetdomains = NULL);
                                                                     // register a mesh with the given name in the bucket

    Mesh_ptr fetch_mesh(const std::string &name);                    // return a (std shared) pointer to a mesh with the given name
    
    MeshFunction_size_t_ptr fetch_celldomains(const std::string &name);     // return a (std shared) pointer to a meshfunction with the given name

    MeshFunction_size_t_ptr fetch_facetdomains(const std::string &name);    // return a (std shared) pointer to a meshfunction with the given name
    
    Mesh_it meshes_begin();                                          // return an iterator to the beginning of the meshes

    Mesh_const_it meshes_begin() const;                              // return a constant iterator to the beginning of the meshes

    Mesh_it meshes_end();                                            // return an iterator to the end of the meshes

    Mesh_const_it meshes_end() const;                                // return a constant iterator to the end of the meshes
 
    MeshFunction_size_t_it celldomains_begin();                      // return an iterator to the beginning of the meshfunctions

    MeshFunction_size_t_const_it celldomains_begin() const;          // return a constant iterator to the beginning of the meshfunctions

    MeshFunction_size_t_it celldomains_end();                        // return an iterator to the end of the meshefunctions

    MeshFunction_size_t_const_it celldomains_end() const;            // return a constant iterator to the end of the meshfunctions
 
    MeshFunction_size_t_it facetdomains_begin();                      // return an iterator to the beginning of the meshfunctions

    MeshFunction_size_t_const_it facetdomains_begin() const;          // return a constant iterator to the beginning of the meshfunctions

    MeshFunction_size_t_it facetdomains_end();                        // return an iterator to the end of the meshefunctions

    MeshFunction_size_t_const_it facetdomains_end() const;            // return a constant iterator to the end of the meshfunctions
 
    void register_visfunctionspace(
               FunctionSpace_ptr functionspace, Mesh_ptr mesh);      // register a visualization functionspace with a mesh

    FunctionSpace_ptr fetch_visfunctionspace(const Mesh_ptr mesh);   // return a (std shared) pointer to a functionspace on the given mesh
    
    Mesh_FunctionSpace_it visfunctionspaces_begin();                 // return an iterator to the beginning of the visfunctionspaces

    Mesh_FunctionSpace_const_it visfunctionspaces_begin() const;     // return a constant iterator to the beginning of the visfunctionspaces

    Mesh_FunctionSpace_it visfunctionspaces_end();                   // return an iterator to the end of the visfunctionspaces

    Mesh_FunctionSpace_const_it visfunctionspaces_end() const;       // return a constant iterator to the end of the visfunctionspaces
 
    //***************************************************************|***********************************************************//
    // System data access
    //***************************************************************|***********************************************************//

    void register_system(SystemBucket_ptr system,                    // register a system with the given name in the bucket
                                          const std::string &name);

    SystemBucket_ptr fetch_system(const std::string &name);          // return a (std shared) pointer to a system with the given name
    
    const SystemBucket_ptr fetch_system(const std::string &name)     // return a constant (std shared) pointer to a system with the given name 
          const;                                                
    
    SystemBucket_it systems_begin();                                 // return an iterator to the beginning of the systems

    SystemBucket_const_it systems_begin() const;                     // return a constant iterator to the beginning of the systems

    SystemBucket_it systems_end();                                   // return an iterator to the end of the systems

    SystemBucket_const_it systems_end() const;                       // return a constant iterator to the end of the systems

    //***************************************************************|***********************************************************//
    // UFL symbol data access
    //***************************************************************|***********************************************************//

    void register_baseuflsymbol(const std::string &baseuflsymbol,    // register a function name with a ufl symbol
                                const std::string &uflsymbol);  

    const std::string fetch_baseuflsymbol(const std::string          // return a string containing the function name belonging to a ufl symbol
                                                &uflsymbol) const;

    const bool contains_baseuflsymbol(const std::string &uflsymbol)  // return a boolean, true if the bucket contains a base ufl symbol
                                                        const;       //                 for the given ufl symbol false otherwise

    void register_uflsymbol(const std::pair< std::string,            // register a (std shared) pointer to a function with a ufl  
                            GenericFunction_ptr > &uflfunctionpair); // symbol

    void register_uflsymbol(GenericFunction_ptr function,            // register a (std shared) pointer to a function with a ufl
                            const std::string &uflsymbol);           // symbol

    GenericFunction_ptr fetch_uflsymbol(const std::string &uflsymbol)// fetch a (std shared) pointer to a function associated with
                       const;                                        // a ufl symbol

    //***************************************************************|***********************************************************//
    // Coefficient functionspace data access
    //***************************************************************|***********************************************************//

    void register_coefficientspace(FunctionSpace_ptr                 // register a functionspace for a coefficient function with a 
                           coefficientspace,                         // given ufl symbol
                           const std::string &uflsymbol);

    const bool contains_coefficientspace(const std::string           // returns a boolean, true if the bucket contains a
                                               &uflsymbol) const;    // functionspace for the uflsymbol, false otherwise

    FunctionSpace_ptr fetch_coefficientspace(const std::string       // returns a (std shared) pointer to a functionspace for
                                               &uflsymbol) const;    // the coefficient with the given uflsymbol

    //***************************************************************|***********************************************************//
    // Detector data access
    //***************************************************************|***********************************************************//

    void register_detector(GenericDetectors_ptr detector,            // register a detector set with the given name
                                  const std::string &name);
    
    GenericDetectors_ptr fetch_detector(const std::string &name);    // fetch a detector set with the given name

    GenericDetectors_it detectors_begin();                           // return an iterator to the beginning of the detectors
    
    GenericDetectors_const_it detectors_begin() const;               // return a constant iterator to the beginning of the detectors
    
    GenericDetectors_it detectors_end();                             // return an iterator to the beginning of the detectors
    
    GenericDetectors_const_it detectors_end() const;                 // return a constant iterator to the beginning of the detectors

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    XDMFFile_ptr fetch_visfile(const Mesh_ptr mesh, bool &newfile);  // fetch (and possibly allocate) a visualization file

    XDMFFile_ptr fetch_convvisfile(const Mesh_ptr mesh, bool &newfile);// fetch (and possibly allocate) a visualization file

    const bool write_vischeckpoints() const                          // return if we're visualizing using checkpoint format or not
    { return write_vischeckpoints_; }

    void output(const int &location);                                // output diagnostics for the bucket

    void checkpoint(const int &location);                            // work out if we're checkpointing the bucket

    const std::string str() const;                                   // return a string describing the bucket contents
    
    virtual const std::string meshes_str(const int &indent=0) const; // return an indented string describing the meshes in the bucket

    const std::string systems_str(const int &indent=0) const;        // return an indented string describing the systems in the bucket

    virtual const std::string coefficientspaces_str(const int        // return an indented string describing the coefficient functionspaces
                                                    &indent=0) const;// contained in the bucket

    virtual const std::string uflsymbols_str(const int &indent=0)    // return an indented string describing the ufl symbols contained in the
                                                           const;    // bucket

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    int dimension_;                                                  // geometry dimension
                                                                     // (assumed size of various objects that don't state
                                                                     //  their own size explicitly)
    
    double_ptr start_time_, old_time_,
                            current_time_, finish_time_, 
                            walltime_limit_;                         // the current and finish times of the simulation

    int_ptr timestep_count_, number_timesteps_;                      // the number of timesteps and number of nonlinear iterations taken

    std::pair< std::string, Constant_ptr > timestep_;                // the timestep, represented as a dolfin constant so it can be used in
                                                                     // the ufl (ufl symbol first member of pair)

    double_ptr timestep_increasetol_;                                // increase tolerance for the timestep

    std::vector< std::pair< FunctionBucket_ptr, double > >           // constraints on the timestep in < FunctionBucket_ptr, max value > pairs
                                              timestep_constraints_;

    int_ptr iteration_count_;                                        // the number of iterations requested and the number of nonlinear 
                                                                     // iterations taken
    int minits_, maxits_;                                            // nonlinear system iteration counts

    double *rtol_, atol_;                                            // nonlinear system tolerances

    double relax_;                                                   // relaxation parameter

    bool ignore_failures_;                                           // ignore convergence failures of the nonlinear systems

    double_ptr steadystate_tol_;                                     // the steady state tolerance

    std::string output_basename_;                                    // the output base name

    bool write_vischeckpoints_;                                      // whether to write visualization checkpoints or not
    
    double_ptr visualization_period_, statistics_period_,            // dump periods
               steadystate_period_, detectors_period_,
               timestepadapt_period_, checkpoint_period_;

    double_ptr visualization_dumptime_, statistics_dumptime_,        // dump periods
               steadystate_dumptime_, detectors_dumptime_,
               timestepadapt_time_, checkpoint_time_;

    int_ptr visualization_period_timesteps_,                         // dump periods in timesteps
            statistics_period_timesteps_,
            steadystate_period_timesteps_, 
            detectors_period_timesteps_,
            timestepadapt_period_timesteps_,
            checkpoint_period_timesteps_;

    int_ptr checkpoint_count_;                                       // the checkpoint count

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    ordered_map<const std::string, Mesh_ptr> meshes_;                // a map from mesh names to (std shared) pointers to meshes

    ordered_map<const std::string, MeshFunction_size_t_ptr> celldomains_;   // a map from mesh names to (std shared) pointers to
                                                                     // meshfunctions describing the cellids

    ordered_map<const std::string, MeshFunction_size_t_ptr> facetdomains_;  // a map from mesh names to (std shared) pointers to
                                                                     // meshfunctions describing the facetids

    std::map< Mesh_ptr, FunctionSpace_ptr > visfunctionspaces_;      // pointers to visualization functionspaces

    ordered_map<const std::string, SystemBucket_ptr> systems_;             // a map from system names to (std shared) pointers to systems

    ordered_map<const std::string, GenericDetectors_ptr> detectors_;        // a map from detector set name to (std shared) pointers to detectors

    //***************************************************************|***********************************************************//
    // Diagnostics data
    //***************************************************************|***********************************************************//

    StatisticsFile_ptr statfile_;                                    // pointer to a statistics file

    DetectorsFile_ptr detfile_;                                      // pointer to a detectors file

    SteadyStateFile_ptr steadyfile_;                                 // pointer to a steady state file

    SystemsConvergenceFile_ptr convfile_;                            // nonlinear systems convergence file

    std::map< Mesh_ptr, XDMFFile_ptr > visfiles_, convvisfiles_;     // pointer to visualization file(s)

    bool write_convvis_;                                             // write convvisfiles_ every nonlinear systems iteration

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill_uflsymbols_();                                         // fill the ufl symbol data structures

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to members of this class

    //***************************************************************|***********************************************************//
    // Base data (continued)
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the name of the bucket

    static time_t start_walltime_;                                   // the start time                                    

    static boost::timer::cpu_timer timer_;                           // timer from the start of the simulation (init)

    //***************************************************************|***********************************************************//
    // Pointers data (continued)
    //***************************************************************|***********************************************************//

    std::map< std::string, std::string > baseuflsymbols_;   // a map from derived ufl symbols to their base symbol
    
    std::map< std::string, GenericFunction_ptr > uflsymbols_;        // a map from derived ufl symbols to field and coefficient pointers

    std::map< std::string, FunctionSpace_ptr > coefficientspaces_;   // a map from the base ufl symbol to a coefficient
                                                                     // functionspace

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    bool steadystate_();                                             // check if a steady state has been attained

    bool perform_action_(double_ptr action_period, 
                         double_ptr previous_action_time, 
                         int_ptr    action_period_timestep,
                         bool       default_action=true);            // indicate if an action should be performed based on periods

    void solve_at_start_();                                          // solve the solvers in this system (in order at the start of a
                                                                     // simulation)

    void solve_in_timeloop_();                                       // solve the solvers in this system (in order during the
                                                                     // timeloop of a simulation)

    bool complete_iterating_(const double &aerror0);                 // indicate if nonlinear systems iterations are complete or not

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    void checkpoint_(const double_ptr time);                         // checkpoint the bucket

    virtual void checkpoint_options_(const double_ptr time);         // checkpoint the options system for the bucket

  };

  typedef std::shared_ptr< Bucket > Bucket_ptr;                    // define a std shared ptr type for the class

}
#endif
