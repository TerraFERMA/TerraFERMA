
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
    
    // fill in the information related to the functionals of this function
    void functionals_fill_(const std::string &optionpath);

  public:
    
    // Specific constructor - with no iterated or old functions
    SpudFunctionBucket(std::string name, std::string optionpath, 
                       GenericFunction_ptr function, System* system);
    
    // Specific constructor
    SpudFunctionBucket(std::string name, std::string optionpath, 
                       GenericFunction_ptr function, GenericFunction_ptr oldfunction,
                       GenericFunction_ptr iteratedfunction, System* system);
    
    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudFunctionBucket();

    // Fill the function assuming the buckettools schema
    void fill();

    // Register a functional in the function (with a spud optionpath)
    void register_functional(Form_ptr functional, std::string name, std::string optionpath);

    // Return a string object describing the function
    std::string str() const
    { str(0); }

    // Return a string object describing the function
    std::string str(int indent) const;

    // Return a string object describing the functionals
    std::string functionals_str() const
    { str(0); }

    // Return a string object describing the functionals
    std::string functionals_str(int indent) const;

    // Return the base optionpath for this function
    std::string optionpath() const
    { return optionpath_; }
    
  };
 
  // Define a boost shared pointer for a spud function
  typedef boost::shared_ptr< SpudFunctionBucket > SpudFunctionBucket_ptr;

}
#endif
