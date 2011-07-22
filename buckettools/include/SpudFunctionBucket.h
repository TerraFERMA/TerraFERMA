
#ifndef __SPUD_FUNCTIONBUCKET_H
#define __SPUD_FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The SpudFunctionBucket class is a derived class of the function that populates the
  // data structures within a function using the spud options parser (and assumes
  // the structure of the buckettools schema)
  class SpudFunctionBucket : public FunctionBucket
  {
  // only accessible by this class
  private:

    // supplement the base class with an optionpath
    std::string optionpath_;

    // supplementary to the function base data store the optionpaths for the functionals_
    std::map< std::string, std::string > functional_optionpaths_;
    
    // fill in the base data
    void base_fill_(const uint &index);
 
    // fill in the base data that aliased coefficients don't have
    void nonaliased_base_fill_(const std::string &optionpath);
 
    void initialize_field_();

    void ic_fill_(const std::string &optionpath);

    void bc_component_fill_(const std::string &optionpath,
                            const std::string &bcname,
                            const std::vector<int> &bcids,
                            const MeshFunction_uint_ptr &edgeidmeshfunction);

    void bc_fill_(const std::string &optionpath,
                  const MeshFunction_uint_ptr &edgeidmeshfunction);

    void initialize_expression_coeff_();

    // fill in the information related to the functionals of this function
    void functionals_fill_();

  public:
    
    // Specific constructor
    SpudFunctionBucket(std::string optionpath, SystemBucket* system);
    
    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudFunctionBucket();

    // Fill the function assuming the buckettools schema
    void field_fill(const uint &index);

    // Fill the function assuming the buckettools schema
    void coeff_fill(const uint &index);

    void initialize_function_coeff();

    void initialize_aliased_coeff();

    // Register a functional in the function (with a spud optionpath)
    void register_functional(Form_ptr functional, std::string name, std::string optionpath);

    // Return a string object describing the function
    const std::string str() const
    { return str(0); }

    // Return a string object describing the function
    const std::string str(int indent) const;

    // Return a string object describing the functionals
    const std::string functionals_str() const
    { return str(0); }

    // Return a string object describing the functionals
    const std::string functionals_str(int indent) const;

    // Return the base optionpath for this function
    const std::string optionpath() const
    { return optionpath_; }

    const bool include_in_diagnostics() const;
    
  };
 
  // Define a boost shared pointer for a spud function
  typedef boost::shared_ptr< SpudFunctionBucket > SpudFunctionBucket_ptr;

}
#endif
