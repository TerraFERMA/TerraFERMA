
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>

namespace buckettools
{

  class Bucket;                                                      // predeclaration
  typedef boost::shared_ptr< Bucket > Bucket_ptr;                    // so we can predeclare a pointer to it
  class SystemBucket;                                                // predeclare the class itself
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  
  enum solve_location { SOLVE_START, SOLVE_TIMELOOP, SOLVE_DIAGNOSTICS };

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

    const double maxchange();                                        // return the maximum change in the (requested) fields over the last timestep

    void updatechange();                                             // update the change function

    void resetchange();                                              // reset the change boolean

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_and_initialize();                                    // attach functions to the forms and functionals
                                                                     // in the system and initialize the matrices

    virtual void copy_diagnostics(SystemBucket_ptr &system, 
                            Bucket_ptr &bucket) const;               // copy the data necessary for the diagnostics data file(s)

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return the name of this system
    { return name_; }

    const std::string uflsymbol() const                              // return the system ufl symbol
    { return uflsymbol_; }

    const Mesh_ptr mesh() const                                      // return a (boost shared) pointer to the system mesh
    { return mesh_; }

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

    Bucket* bucket()                                                 // return a pointer to the parent bucket
    { return bucket_; }

    const Bucket* bucket() const                                     // return a constant pointer to the parent bucket
    { return bucket_; }

    const int solve_location() const                                 // return an integer describing where this system is solved
    { return solve_location_; }

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

    //***************************************************************|***********************************************************//
    // Coefficient data access
    //***************************************************************|***********************************************************//

    void register_coeff(FunctionBucket_ptr coeff, 
                                           const std::string &name); // register a coefficient with the given name

    FunctionBucket_ptr fetch_coeff(const std::string &name);         // return a (boost shared) pointer to a coefficient with the
                                                                     // given name

    FunctionBucket_it coeffs_begin();                                // return an iterator to the beginning of the coefficients

    FunctionBucket_const_it coeffs_begin() const;                    // return a constant iterator to the beginning of the
                                                                     // coefficients

    FunctionBucket_it coeffs_end();                                  // return an iterator to the end of the coefficients

    FunctionBucket_const_it coeffs_end() const;                      // return a constant iterator to the end of the coefficients

    //***************************************************************|***********************************************************//
    // Solver bucket data access
    //***************************************************************|***********************************************************//

    void register_solver(SolverBucket_ptr solver, 
                                           const std::string &name); // register a solver bucket with the given name

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

    std::vector<BoundaryCondition_ptr>::iterator bcs_begin();        // return an iterator to the beginning of the system bcs

    std::vector<BoundaryCondition_ptr>::const_iterator bcs_begin()   // return a constant iterator to the beginning of the system
                                                          const;     // bcs

    std::vector<BoundaryCondition_ptr>::iterator bcs_end();          // return an iterator to the end of the system bcs

    std::vector<BoundaryCondition_ptr>::const_iterator bcs_end()     // return a constant iterator to the end of the system bcs
                                                          const;

    const std::vector< BoundaryCondition_ptr > bcs() const           // return a constant vector of system bcs
    { return bcs_; }
    
    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    void output(const bool &write_vis);                              // output the diagnostics on this system

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

    virtual const std::string fields_str(const int &indent) const;          // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string coeffs_str() const                     // return a string describing the coefficients in the system
    { return coeffs_str(0); }

    virtual const std::string coeffs_str(const int &indent) const;          // return an indented string describing the fields in the
                                                                     // system

    virtual const std::string solvers_str() const                    // return a string describing the solver buckets in the system
    { return coeffs_str(0); }

    virtual const std::string solvers_str(const int &indent) const;         // return an indented string describing the solver buckets in
                                                                     // the system

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

    void collect_bcs_();                                             // collect a vector of (boost shared) pointers bcs from the
                                                                     // fields

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the system name

    std::string uflsymbol_;                                          // the system ufl symbol

    Bucket* bucket_;                                                 // a pointer to the parent bucket

    Mesh_ptr mesh_;                                                  // a (boost shared) pointer to the system mesh

    FunctionSpace_ptr functionspace_;                                // a (boost shared) pointer to the system functionspace

    Function_ptr function_, oldfunction_, iteratedfunction_;         // (boost shared) pointers to the system functions at different
                                                                     // time levels (old, iterated - most up to date -, base)

    int solve_location_;                                             // when this system will be solved

    Function_ptr changefunction_;                                    // (boost shared) pointer to the change between timesteps

    bool_ptr change_calculated_;                                     // indicate if the change has been recalculated recently

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, FunctionBucket_ptr > fields_;             // a map from field names to (boost shared) pointers to fields
    
    std::map< std::string, FunctionBucket_ptr > coeffs_;             // a map from coefficient names to (boost shared) pointers to
                                                                     // coefficients
    
    std::map< std::string, SolverBucket_ptr > solvers_;              // a map from solver bucket names to (boost shared) pointers to
                                                                     // solver buckets

    std::map< int, SolverBucket_ptr > orderedsolvers_;               // an ordered (user defined) map from solver bucket names to
                                                                     // (boost shared) pointers to solver buckets

    std::vector< BoundaryCondition_ptr > bcs_;                       // a vector of (boost shared) poitners to bcs

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the data structures in this system

  };

  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // define a (boost shared) pointer to the system class type

}
#endif
