
#ifndef __SOLVERBUCKET_H
#define __SOLVERBUCKET_H

#include "BoostTypes.h"
#include "BucketPETScBase.h"
#include <dolfin.h>
#include "petscsnes.h"

namespace buckettools
{

  class SystemBucket;                                                // predeclare
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  class FunctionBucket;                                              // predeclare the class itself
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // so we can predeclare a pointer to it
  
  //*****************************************************************|************************************************************//
  // SolverBucket class:
  //
  // The SolverBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  //*****************************************************************|************************************************************//
  class SolverBucket
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

    void solve();                                                    // run the nonlinear solver described by this class

    void assemble_linearforms(const bool &reset_tensor);             // assemble all linear forms in this solver

    void assemble_bilinearforms(const bool &reset_tensor);           // assemble all bilinear forms in this solver

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_form_coeffs();                                       // attach coefficients to the forms in this solver

    void initialize_matrices();                                      // initialize the vectors and matrices (preassembly) described
                                                                     // by the forms in this solver

    virtual void copy_diagnostics(SolverBucket_ptr &solver, 
                                  SystemBucket_ptr &system) const;   // copy the data necessary for the diagnostics data file(s)

    void initialize_diagnostics() const;                             // initialize any diagnostic output in the solver

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a string containing the solver name
    { return name_; }

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

    const bool kspnullspace_monitor() const;                         // return true if we're using a visualization monitor
    
    const ConvergenceFile_ptr convergence_file() const;              // return a pointer to the convergence file

    const KSPConvergenceFile_ptr ksp_convergence_file() const;       // return a pointer to the ksp convergence file

    const bool monitor_norms() const                                 // return true if norms should be monitored in nonlinear iterations
    { return monitornorms_; }

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
    // Output functions
    //***************************************************************|***********************************************************//

    virtual const std::string str() const                            // return a string describing the contents of this solver
    { return str(0); }

    virtual const std::string str(int indent) const;                 // return an indented string describing the contents of this
                                                                     // solver

    virtual const std::string forms_str() const                      // return a string describing the forms in this solver
    { return forms_str(0); }

    virtual const std::string forms_str(const int &indent) const;    // return an indented string describing the forms in this
                                                                     // solver

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

    PETScMatrix_ptr matrix_, matrixpc_;                              // dolfin petsc matrix types

    PETScVector_ptr rhs_, res_, work_;                               // dolfin petsc vector types

    double rtol_, atol_, stol_;                                      // nonlinear solver tolerances

    int minits_, maxits_, maxfes_;                                   // nonlinear solver iteration counts

    int_ptr iteration_count_;                                        // nonlinear iterations taken (may not be accurate!)

    bool ident_zeros_, ident_zeros_pc_;                              // replace zero rows with the identity (matrix and pc)

    bool ignore_failures_;                                           // ignore solver failures

    std::string name_;                                               // solver name

    std::string type_;                                               // solver type (string)

    SystemBucket* system_;                                           // parent system

    CustomMonitorCtx snesmctx_, kspmctx_;                            // monitor contexts

    ConvergenceFile_ptr convfile_;                                   // diagnostic convergence file

    KSPConvergenceFile_ptr kspconvfile_;                             // diagnostic convergence file

    bool_ptr kspnullspacemonitor_;                                   // visualization monitors

    bool copy_;                                                      // flag if this is a diagnostic copy or not

    bool monitornorms_;                                              // monitor the norms in nonlinear iterations

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, Form_ptr > forms_;                        // a map from the form names to the form pointers

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the class data structures

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

  typedef boost::shared_ptr< SolverBucket > SolverBucket_ptr;        // define a (boost shared) pointer to the solver class type

}
#endif
