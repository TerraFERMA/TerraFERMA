
#ifndef __FUNCTIONBUCKET_H
#define __FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  class SystemBucket;                                                // predeclaration
  
  //*****************************************************************|************************************************************//
  // FunctionBucket class:
  //
  // The FunctionBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  //*****************************************************************|************************************************************//
  class FunctionBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    FunctionBucket();                                                // default constructor

    FunctionBucket(SystemBucket* system);                            // specific constructor
    
    ~FunctionBucket();                                               // default destructor

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a constant string giving the function name
    { return name_; }

    const std::string type() const                                   // return a constant string giving the function type
    { return type_; }

    const std::string uflsymbol() const                              // return a constant string giving the ufl symbol 
    { return uflsymbol_; }                                           // for this function

    const uint index() const                                         // return a constant unsigned integer to the index of this
    { return index_; }                                               // function in the parent system

    const SystemBucket* system() const                               // return a constant pointer to the parent system
    { return system_; }

    const FunctionSpace_ptr functionspace() const                    // return a constant (boost shared) pointer to the
    { return functionspace_; }                                       // functionspace
                                                                     // NOTE: if this is a field of a mixed system functionspace,
                                                                     // this will return a subspace


    const GenericFunction_ptr function() const                       // return a constant (boost shared) pointer to the 
    { return function_; }                                            // function 
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr oldfunction() const                    // return a constant (boost shared) pointer to the old
    { return oldfunction_; }                                         // function (previous timestep's values)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    // Return the iteratedfunction_
    const GenericFunction_ptr iteratedfunction() const               // return a constant (boost shared) pointer to the iterated
    { return iteratedfunction_; }                                    // function (most up to date values within an iteration)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    // Return the icexpression_
    const Expression_ptr icexpression() const                        // return a constant (boost shared) pointer to the initial
    { return icexpression_; }                                        // condition expression for this function

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_functional_coeffs();                                 // attach the coefficients to the functionals of this function

    //***************************************************************|***********************************************************//
    // Functional data access
    //***************************************************************|***********************************************************//

    void register_functional(Form_ptr functional, std::string name); // register a functional with the given name in this function

    Form_ptr fetch_functional(std::string name);                     // return a (boost shared) pointer to a functional with the given name

    Form_it functionals_begin();                                     // return an iterator to the beginning of the functionals of this function

    Form_const_it functionals_begin() const;                         // return a constant iterator to the beginning of the functionals of this function

    Form_it functionals_end();                                       // return an iterator to the end of the functionals of this function

    Form_const_it functionals_end() const;                           // return a constant iterator to the end of the functionals of this function

    //***************************************************************|***********************************************************//
    // BC data access
    //***************************************************************|***********************************************************//

    void register_bcexpression(Expression_ptr bcexpression,          // register an expression for a bc in this function
                                               std::string name);

    void register_bc(BoundaryCondition_ptr bc, std::string name);    // register a bc in this function

    BoundaryCondition_it bcs_begin();                                // return an iterator to the beginning of the bcs of this function

    BoundaryCondition_const_it bcs_begin() const;                    // return a constant iterator to the beginning of the bcs of this function

    BoundaryCondition_it bcs_end();                                  // return an iterator to the end of the bcs of this function

    BoundaryCondition_const_it bcs_end() const;                      // return a constant iterator to the end of the bcs of this function

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    void output();                                                   // output diagnostics about this function

    virtual const bool include_in_diagnostics() const;               // return a boolean indicating if this function is included in 
                                                                     // diagnostic output

    virtual const std::string str() const                            // return a string describing the contents of this function
    { return str(0); }

    virtual const std::string str(int indent) const;                 // return an indented string describing the contents 
                                                                     // of this function

    virtual const std::string functionals_str() const                // return a string describing the functionals of this function
    { return functionals_str(0); }

    virtual const std::string functionals_str(int indent) const;     // return an indented string describing the functionals 
                                                                     // of this function

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the data structures of this function bucket

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the function name

    std::string uflsymbol_;                                          // the function ufl symbol

    uint index_;                                                     // the index of this field or coefficient in the system
                                                                     // (most relevant for field which are subfunctions of the system function)

    SystemBucket* system_;                                           // a pointer to the parent system

    FunctionSpace_ptr functionspace_;                                // the functionspace (may be a subspace) of this function (if it is
                                                                     // a field or a coefficient function)

    GenericFunction_ptr function_, oldfunction_, iteratedfunction_;  // (boost shared) pointers to different timelevel values of this function

    Expression_ptr icexpression_;                                    // (boost shared) pointer to an expression describing the initial condition

    int size_;                                                       // size of the function (most relevant for rank 1, vectors)

    std::vector< int > shape_;                                       // shape of the function (most relevant for rank 2, tensors)

    std::string rank_;                                               // a *string* describing the rank of the function
    
    std::string type_;                                               // a *string* describing the type of function (function, expression, constant)
    
    File_ptr pvdfile_;                                               // (boost shared) pointer to a pvd file output for this function

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, Form_ptr > functionals_;                  // map from functional names to form (boost shared) pointers

    std::map< std::string, Expression_ptr > bcexpressions_;          // map from bc names to bc expression (boost shared) pointers
    
    std::map< std::string, BoundaryCondition_ptr > bcs_;             // map from bc::id names to (boost shared) pointers to bcs
    
  };

  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // define a (boost shared) pointer to the function bucket class type

}
#endif
