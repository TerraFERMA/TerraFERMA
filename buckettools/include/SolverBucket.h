
#ifndef __SOLVERBUCKET_H
#define __SOLVERBUCKET_H

#include "BoostTypes.h"
#include "BucketPETScBase.h"
#include <dolfin.h>
#include "petscsnes.h"

namespace buckettools
{

  class SystemBucket;                                                // predeclaration
  
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

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a string containing the solver name
    { return name_; }

    const std::string type() const                                   // return a string describing the solver type
    { return type_; }

    //***************************************************************|***********************************************************//
    // Form data access
    //***************************************************************|***********************************************************//

    void register_form(Form_ptr form, std::string name);             // register a form in the solver

    bool contains_form(std::string name);                            // return a boolean, true if the names form exists in the
                                                                     // solver, false otherwise

    Form_ptr fetch_form(std::string name);                           // fetch the named form

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

    virtual const std::string forms_str(int indent) const;           // return an indented string describing the forms in this
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

    std::string name_;                                               // solver name

    std::string type_;                                               // solver type (string)

    SystemBucket* system_;                                           // parent system

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, Form_ptr > forms_;                        // a map from the form names to the form pointers

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the class data structures

  };

  typedef boost::shared_ptr< SolverBucket > SolverBucket_ptr;        // define a (boost shared) pointer to the solver class type

}
#endif
