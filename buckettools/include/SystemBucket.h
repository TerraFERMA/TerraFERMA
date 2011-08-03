
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>

namespace buckettools
{

  // Predeclaration of parent classes to allow two-way dependencies
  class Bucket;
  
  // The SystemBucket class describes a functionspace and a set of solvers that act on the fields
  // contained in that (potentially mixed) functionspace.
  // This base class describes the basic data structures while derived classes may be defined
  // that allow it to be linked to an options system.
  class SystemBucket
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
    
    std::map< std::string, SolverBucket_ptr > solvers_;

    std::map< int, SolverBucket_ptr > orderedsolvers_;

    std::vector< BoundaryCondition_ptr > bcs_;

    void attach_all_coeffs_();

    void attach_function_coeffs_(FunctionBucket_it f_begin, FunctionBucket_it f_end);

    void attach_solver_coeffs_(SolverBucket_it s_begin, SolverBucket_it s_end);

    void collect_bcs_();

  public:

    // Default constructor
    SystemBucket();

    // Specific constructor
    SystemBucket(Bucket* bucket);
    
    // Default destructor
    ~SystemBucket();

    // Register a field (subfunction) in the system
    void register_field(FunctionBucket_ptr field, std::string name);

    // Fetch a field (subfunction) from the system
    FunctionBucket_ptr fetch_field(std::string name);

    const FunctionBucket_ptr fetch_field(std::string name) const;

    FunctionBucket_it fields_begin();

    FunctionBucket_const_it fields_begin() const;

    FunctionBucket_it fields_end();

    FunctionBucket_const_it fields_end() const;

    // Register a coefficient (expression or function) in the system
    void register_coeff(FunctionBucket_ptr coeff, std::string name);

    // Fetch a coefficient (expression or function) from the system
    FunctionBucket_ptr fetch_coeff(std::string name);

    FunctionBucket_it coeffs_begin();

    FunctionBucket_const_it coeffs_begin() const;

    FunctionBucket_it coeffs_end();

    FunctionBucket_const_it coeffs_end() const;

    void register_solver(SolverBucket_ptr solver, std::string name);

    SolverBucket_it solvers_begin();

    SolverBucket_const_it solvers_begin() const;

    SolverBucket_it solvers_end();

    SolverBucket_const_it solvers_end() const;

    int_SolverBucket_it orderedsolvers_begin();

    int_SolverBucket_const_it orderedsolvers_begin() const;

    int_SolverBucket_it orderedsolvers_end();

    int_SolverBucket_const_it orderedsolvers_end() const;

    std::vector<BoundaryCondition_ptr>::iterator bcs_begin();

    std::vector<BoundaryCondition_ptr>::const_iterator bcs_begin() const;

    std::vector<BoundaryCondition_ptr>::iterator bcs_end();

    std::vector<BoundaryCondition_ptr>::const_iterator bcs_end() const;

    const std::vector< BoundaryCondition_ptr > bcs() const
    { return bcs_; }
    
    // Return a string describing the contents of the system
    virtual const std::string str() const
    { return str(0); }

    // Return a string describing the contents of the system
    virtual const std::string str(int indent) const;

    // Print a description of the fields contained in the system
    virtual const std::string fields_str() const
    { return fields_str(0); }

    // Print a description of the fields contained in the system
    virtual const std::string fields_str(int indent) const;

    // Print a description of the coefficients contained in the system
    virtual const std::string coeffs_str() const
    { return coeffs_str(0); }

    // Print a description of the coefficients contained in the system
    virtual const std::string coeffs_str(int indent) const;

    // Print a description of the solvers contained in the system
    virtual const std::string solvers_str() const
    { return coeffs_str(0); }

    // Print a description of the solvers contained in the system
    virtual const std::string solvers_str(int indent) const;

    void attach_and_initialize();

    void solve();

    // Return the system name
    const std::string name() const
    { return name_; }

    // Return the system uflsymbol
    const std::string uflsymbol() const
    { return uflsymbol_; }

    // Return a pointer to the system mesh
    const Mesh_ptr mesh() const
    { return mesh_; }

    const FunctionSpace_ptr functionspace() const
    { return functionspace_; }

    const Function_ptr function() const
    { return function_; }

    const Function_ptr oldfunction() const
    { return oldfunction_; }

    const Function_ptr iteratedfunction() const
    { return iteratedfunction_; }

    Bucket* bucket()
    { return bucket_; }

    const Bucket* bucket() const
    { return bucket_; }

    void output();

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;

}
#endif
