
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>

namespace buckettools
{

  // Predeclaration of parent classes to allow two-way dependencies
  class Bucket;
  
  // The System class describes a functionspace and a set of solvers that act on the fields
  // contained in that (potentially mixed) functionspace.
  // This base class describes the basic data structures while derived classes may be defined
  // that allow it to be linked to an options system.
  class System
  {
  // only accessible by this class
  private:

    // Empty the data structures of the system
    void empty_();

  protected:

    // the system name
    std::string name_;

    // the system uflsymbol
    std::string uflsymbol_;

    // a pointer to the parent bucket of this system
    Bucket* bucket_;

    // a pointer to the mesh which this system is on
    Mesh_ptr mesh_;

    // a pointer to the functionspace which this system describes
    FunctionSpace_ptr functionspace_;

    // the function (at several time-levels) on the above functionspace_
    Function_ptr function_, oldfunction_, iteratedfunction_;

    // a map from field names to the functionbuckets of fields (subfunctions)
    std::map< std::string, FunctionBucket_ptr > fields_;
    
    // a map from coefficient names to the functionbuckets of coefficients
    std::map< std::string, FunctionBucket_ptr > coeffs_;
    
    // a map from ufl symbols to field and coefficient names
    std::map< std::string, std::string > uflnames_;
    
    // a map from ufl symbols to pointers to functions
    std::map< std::string, GenericFunction_ptr > uflsymbols_;

    std::map< std::string, FunctionSpace_ptr > coefficientspaces_;

    std::map< std::string, SolverBucket_ptr > solvers_;
    
  public:

    // Default constructor
    System();

    // Specific constructor
    System(Bucket* bucket);
    
    // Default destructor
    ~System();

    // Register a field (subfunction) in the system
    void register_field(FunctionBucket_ptr field, std::string name);

    // Fetch a field (subfunction) from the system
    FunctionBucket_ptr fetch_field(std::string name);

    // Register a coefficient (expression or function) in the system
    void register_coeff(FunctionBucket_ptr coeff, std::string name);

    // Fetch a coefficient (expression or function) from the system
    FunctionBucket_ptr fetch_coeff(std::string name);

    // Register a ufl name
    void register_uflname(std::string name, std::string uflsymbol);

    // Register a ufl system and function pointer in the system
    void register_uflsymbol(GenericFunction_ptr function, std::string uflsymbol);

    // Return the name of a function with uflsymbol
    std::string fetch_uflname(std::string uflsymbol);

    // Create a ufl system pointing at a null function pointer
    void create_uflsymbol(std::string uflsymbol);

    // Reset a function associated with a ufl system
    void reset_uflsymbol(GenericFunction_ptr function, std::string uflsymbol);

    // Return a pointer to a dolfin GenericFunction with the given uflsymbol
    GenericFunction_ptr fetch_uflsymbol(std::string uflsymbol);

    void register_coefficientspace(FunctionSpace_ptr coefficientspace, std::string name);

    bool contains_coefficientspace(std::string name);

    FunctionSpace_ptr fetch_coefficientspace(std::string name);

    void register_solver(SolverBucket_ptr solver, std::string name);

    // Return a string describing the contents of the system
    virtual std::string str() const
    { str(0); }

    // Return a string describing the contents of the system
    virtual std::string str(int indent) const;

    // Print a description of the fields contained in the system
    virtual std::string uflsymbols_str() const
    { fields_str(0); }

    // Print a description of the fields contained in the system
    virtual std::string uflsymbols_str(int indent) const;

    // Print a description of the fields contained in the system
    virtual std::string fields_str() const
    { fields_str(0); }

    // Print a description of the fields contained in the system
    virtual std::string fields_str(int indent) const;

    // Print a description of the coefficients contained in the system
    virtual std::string coeffs_str() const
    { coeffs_str(0); }

    // Print a description of the coefficients contained in the system
    virtual std::string coeffs_str(int indent) const;

    // Print a description of the solvers contained in the system
    virtual std::string solvers_str() const
    { coeffs_str(0); }

    // Print a description of the solvers contained in the system
    virtual std::string solvers_str(int indent) const;

    // Return the system name
    std::string name() const
    { return name_; }

    // Return the system uflsymbol
    std::string uflsymbol() const
    { return uflsymbol_; }

    // Return a pointer to the system mesh
    Mesh_ptr mesh() const
    { return mesh_; }

    FunctionSpace_ptr functionspace() const
    { return functionspace_; }

    Function_ptr function() const
    { return function_; }

    Function_ptr oldfunction() const
    { return oldfunction_; }

    Function_ptr iteratedfunction() const
    { return iteratedfunction_; }

    Bucket* bucket() const
    { return bucket_; }

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< System > System_ptr;

}
#endif
