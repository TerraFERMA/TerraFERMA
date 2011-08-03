
#ifndef __FUNCTIONBUCKET_H
#define __FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  // Predeclaration of parent classes to allow two-way dependencies
  class SystemBucket;
  
  // The FunctionBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  class FunctionBucket
  {
  // only accessible by this class
  private:

    // Empty the data structures of the function
    void empty_();

  protected:
    
    // the function name
    std::string name_;

    // ufl symbol
    std::string uflsymbol_;

    // the field index
    uint index_;

    // the system to which this function belongs
    SystemBucket* system_;

    FunctionSpace_ptr functionspace_;

    // the function
    GenericFunction_ptr function_, oldfunction_, iteratedfunction_;

    // expression for the initial condition (if it exists)
    Expression_ptr icexpression_;

    File_ptr pvdfile_;

    // a map from functional names to the forms (the real point of this class)
    std::map< std::string, Form_ptr > functionals_;

    // a map from components to the (sub)functionspaces they are described on
    std::map< int, FunctionSpace_ptr > subfunctionspaces_;
    
    // a map from bc names to the Expressions describing them
    std::map< std::string, Expression_ptr > bcexpressions_;
    
    // a map from bc::id names to (mostly dirichlet) boundary conditions
    std::map< std::string, BoundaryCondition_ptr > bcs_;
    
    int size_;

    std::vector< int > shape_;

    // would prefer to switch this to an integer at some point
    std::string rank_;
    
    // would prefer to switch this to an integer at some point
    std::string type_;
    
  public:

    // Default constructor
    FunctionBucket();

    // Specific constructor with an uninitialised name
    FunctionBucket(SystemBucket* system);
    
    // Default destructor
    ~FunctionBucket();

    // Register a subfunctionspace in the system
    void register_subfunctionspace(FunctionSpace_ptr subfunctionspace, int component);

    // Return whether a subfunctionspace with the given name is in the system
    bool contains_subfunctionspace(int component);

    // Return a pointer to a dolfin subfunctionspace with the given name
    FunctionSpace_ptr fetch_subfunctionspace(int component);

    // Register a functional in the function
    void register_functional(Form_ptr functional, std::string name);

    // Return a pointer to a functional with the given name
    Form_ptr fetch_functional(std::string name);

    Form_it functionals_begin();

    Form_const_it functionals_begin() const;

    Form_it functionals_end();

    Form_const_it functionals_end() const;

    // Register a bc expression in the functionbucket
    void register_bcexpression(Expression_ptr bcexpression, std::string name);

    // Register a bc expression in the functionbucket
    void register_bc(BoundaryCondition_ptr bc, std::string name);

    BoundaryCondition_it bcs_begin();

    BoundaryCondition_const_it bcs_begin() const;

    BoundaryCondition_it bcs_end();

    BoundaryCondition_const_it bcs_end() const;

    void attach_functional_coeffs();

    // Return a string describing the contents of the function
    virtual const std::string str() const
    { return str(0); }

    // Return a string describing the contents of the function
    virtual const std::string str(int indent) const;

    // Print a description of the functionals contained in the function
    virtual const std::string functionals_str() const
    { return functionals_str(0); }

    // Print a description of the functionals contained in the system
    virtual const std::string functionals_str(int indent) const;

    // Return the function name
    const std::string name() const
    { return name_; }

    // Return the function type
    const std::string type() const
    { return type_; }

    // Return the function uflsymbol
    const std::string uflsymbol() const
    { return uflsymbol_; }

    // Return the function index
    const uint index() const
    { return index_; }

    const SystemBucket* system() const
    { return system_; }

    // Return the function_
    const FunctionSpace_ptr functionspace() const
    { return functionspace_; }

    // Return the function_
    const GenericFunction_ptr function() const
    { return function_; }

    // Return the oldfunction_
    const GenericFunction_ptr oldfunction() const
    { return oldfunction_; }

    // Return the iteratedfunction_
    const GenericFunction_ptr iteratedfunction() const
    { return iteratedfunction_; }

    // Return the icexpression_
    const Expression_ptr icexpression() const
    { return icexpression_; }

    virtual const bool include_in_diagnostics() const;

    void output();

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;

}
#endif
