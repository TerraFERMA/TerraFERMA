
#ifndef __SOLVERBUCKET_H
#define __SOLVERBUCKET_H

#include "BoostTypes.h"
#include <dolfin.h>
#include "petscsnes.h"

namespace buckettools
{

  // Predeclaration of parent classes to allow two-way dependencies
  class SystemBucket;
  
  // The SolverBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
  class SolverBucket
  {
  // only accessible by this class
  private:

    // Empty the data structures of the function
    void empty_();

  protected:
    
    SNES snes_;

    KSP ksp_;

    // the function name
    std::string name_;

    // the function type
    std::string type_;

    // the field index
    uint index_;

    // the system to which this function belongs
    SystemBucket* system_;

    // a map from form names to the forms
    std::map< std::string, Form_ptr > forms_;

  public:

    // Default constructor
    SolverBucket();

    // Specific constructor with an uninitialised name
    SolverBucket(SystemBucket* system);
    
    // Default destructor
    ~SolverBucket();

    // Register a form in the solver
    void register_form(Form_ptr form, std::string name);

    // Return a pointer to a form with the given name
    Form_ptr fetch_form(std::string name);

    Form_it forms_begin();

    Form_const_it forms_begin() const;

    Form_it forms_end();

    Form_const_it forms_end() const;

    void attach_form_coeffs();

    // Return a string describing the contents of the solver
    virtual std::string str() const
    { str(0); }

    // Return a string describing the contents of the solver
    virtual std::string str(int indent) const;

    // Print a description of the forms contained in the solver
    virtual std::string forms_str() const
    { forms_str(0); }

    // Print a description of the forms contained in the solver
    virtual std::string forms_str(int indent) const;

    // Return the function name
    std::string name() const
    { return name_; }

    // Return the function type
    std::string type() const
    { return type_; }

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< SolverBucket > SolverBucket_ptr;

}
#endif
