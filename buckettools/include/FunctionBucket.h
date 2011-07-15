
#ifndef __FUNCTIONBUCKET_H
#define __FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  // Predeclaration of parent classes to allow two-way dependencies
  class System;
  
  // The FunctionBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  class FunctionBucket
  {
  // only accessible by this class
  private:

    // the function name
    std::string name_;

    // Empty the data structures of the function
    void empty_();

  protected:
    
    // the system to which this function belongs
    System* system_;

    // the function
    GenericFunction_ptr function_, oldfunction_, iteratedfunction_;

    // a map from functional names to the forms (the real point of this class)
    std::map< std::string, Form_ptr > functionals_;
    
  public:

    // No default constructor - always require a function pointer

    // Specific constructor with an uninitialised name
    FunctionBucket(GenericFunction_ptr function, GenericFunction_ptr oldfunction, GenericFunction_ptr iteratedfunction, System* system)
    { FunctionBucket("uninitialised_name", function, oldfunction, iteratedfunction, system); }

    // Specific constructor
    FunctionBucket(std::string name, GenericFunction_ptr function, GenericFunction_ptr oldfunction, GenericFunction_ptr iteratedfunction, System* system);
    
    // Default destructor
    ~FunctionBucket();

    // Register a functional in the function
    void register_functional(Form_ptr functional, std::string name);

    // Return a pointer to a functional with the given name
    Form_ptr fetch_functional(std::string name);

    // Return a string describing the contents of the function
    virtual std::string str() const
    { str(0); }

    // Return a string describing the contents of the function
    virtual std::string str(int indent) const;

    // Print a description of the functionals contained in the function
    virtual std::string functionals_str() const
    { functionals_str(0); }

    // Print a description of the functionals contained in the system
    virtual std::string functionals_str(int indent) const;

    // Return the function name
    std::string name() const
    { return name_; }

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;

}
#endif
