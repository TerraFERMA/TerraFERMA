
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

    // the system name
    std::string name_;

    // Empty the data structures of the system
    void empty_();

  protected:

    // a pointer to the parent bucket of this system
    Bucket* bucket_;

    // a pointer to the mesh which this system is on
    Mesh_ptr mesh_;

    // a pointer to the functionspace which this system describes
    FunctionSpace_ptr functionspace_;

    // the function (at several time-levels) on the above functionspace_
    Function_ptr function_, oldfunction_, iteratedfunction_;

    // a map from field names to the subfunctionspaces they are described on
    std::map< std::string, FunctionSpace_ptr > subfunctionspaces_;
    
    // a map from field names to the functionbuckets of fields (subfunctions)
    std::map< std::string, FunctionBucket_ptr > fields_;
    
    // a map from coefficient names to the functionbuckets of coefficients
    std::map< std::string, FunctionBucket_ptr > coeffs_;
    
    // a map from field::bc names to the subfunctionspaces they are described on
    std::map< std::string, Expression_ptr > bcexpressions_;
    
    // a map from field::bc::id names to dirichlet boundary conditions
    std::map< std::string, DirichletBC_ptr > dirichletbcs_;
    
    // a map from component integer index to initial condition expression
    std::map< uint, Expression_ptr > icexpressions_;

    // a map from ufl symbols to field and coefficient names
    std::map< std::string, std::string > uflnames_;
    
    // a map from ufl symbols to pointers to functions
    std::map< std::string, GenericFunction_ptr > uflsymbols_;
    
  public:

    // No default constructor - always require a mesh pointer

    // Specific constructor with an uninitialised name
    System(Mesh_ptr mesh, Bucket* bucket)
    { System("uninitialised_name", mesh, bucket); }

    // Specific constructor
    System(std::string name, Mesh_ptr mesh, Bucket* bucket);
    
    // Default destructor
    ~System();

    // Register a subfunctionspace in the system
    void register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name);

    // Return whether a subfunctionspace with the given name is in the system
    bool contains_subfunctionspace(std::string name);

    // Return a pointer to a dolfin subfunctionspace with the given name
    FunctionSpace_ptr fetch_subfunctionspace(std::string name);

    // Register a field (subfunction) in the system
    void register_field(FunctionBucket_ptr field, std::string name);

    // Register a coefficient (expression or function) in the system
    void register_coeff(FunctionBucket_ptr coeff, std::string name);

    // Register a bc expression in the system
    void register_bcexpression(Expression_ptr bcexpression, std::string name);

    // Register a bc expression in the system
    void register_dirichletbc(DirichletBC_ptr bc, std::string name);

    // Register an initial condition expression
    void register_icexpression(Expression_ptr ic, uint component);

    // Register a ufl name
    void register_uflname(std::string name, std::string uflsymbol);

    // Register a ufl system and function pointer in the system
    void register_uflsymbol(GenericFunction_ptr function, std::string uflsymbol);

    // Create a ufl system pointing at a null function pointer
    void create_uflsymbol(std::string uflsymbol);

    // Reset a function associated with a ufl system
    void reset_uflsymbol(GenericFunction_ptr function, std::string uflsymbol);

    // Return a string describing the contents of the system
    virtual std::string str() const
    { str(0); }

    // Return a string describing the contents of the system
    virtual std::string str(int indent) const;

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

    // Print a description of the bcexpressions contained in the system
    virtual std::string bcexpressions_str() const
    { bcexpressions_str(0); }

    // Print a description of the bcexpressions contained in the system
    virtual std::string bcexpressions_str(int indent) const;

    // Print a description of the icexpressions contained in the system
    virtual std::string icexpressions_str() const
    { icexpressions_str(0); }

    // Print a description of the icexpressions contained in the system
    virtual std::string icexpressions_str(int indent) const;

    // Return the system name
    std::string name() const
    { return name_; }

    // Return a pointer to the system mesh
    Mesh_ptr mesh() const
    { return mesh_; }

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< System > System_ptr;

}
#endif
