
#ifndef __SPUD_FUNCTIONBUCKET_H
#define __SPUD_FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudFunctionBucket class:
  //
  // The SpudFunctionBucket class is a derived class of the function that populates the
  // data structures within a function using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudFunctionBucket : public FunctionBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudFunctionBucket(std::string optionpath, SystemBucket* system);// optional constructor
    
    virtual ~SpudFunctionBucket();                                   // default destructor

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a string containing the optionpath of this function
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void field_fill(const uint &index);                              // fill this function assuming it's a field

    void coeff_fill(const uint &index);                              // fill this function assuming it's a coefficient

    void initialize_function_coeff();                                // initialize this function assuming it's a coefficient function

    //***************************************************************|***********************************************************//
    // Functional data access
    //***************************************************************|***********************************************************//

    void register_functional(Form_ptr functional, std::string name,  // register a functional with the given name and optionpath
                                            std::string optionpath);

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const bool include_in_diagnostics() const;                       // return a boolean indicating if this function is to 
                                                                     // be included in diagnostic output
    
    const std::string str() const                                    // return a string describing the contents of this function
    { return str(0); }

    const std::string str(int indent) const;                         // return an indented string describing the contents of this function

    const std::string functionals_str() const                        // return a string describing the functionals of this function
    { return str(0); }

    const std::string functionals_str(int indent) const;             // return an indented string describing the functionals of this function

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the optionpath of this function

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, std::string > functional_optionpaths_;    // a map from funcional name to functional optionpath
    
    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void base_fill_(const uint &index);                              // fill the base data of this function
 
    void initialize_field_();                                        // initialize this function assuming its a field

    void initialize_expression_coeff_();                             // initialize this function assuming its a coefficient expression

    void ic_fill_(const std::string &optionpath);                    // fill in the initial condition data for this function

    void bc_fill_(const std::string &optionpath,                     // fill in the bc for this function
                  const MeshFunction_uint_ptr 
                                            &edgeidmeshfunction);

    void bc_component_fill_(const std::string &optionpath,           // fill in the bc data for a component of this function
                            const std::string &bcname,
                            const std::vector<int> &bcids,
                            const MeshFunction_uint_ptr 
                                            &edgeidmeshfunction);

    void functionals_fill_();                                        // fill in the data for the functionals of this function

  };
 
  typedef boost::shared_ptr< SpudFunctionBucket >                    // define a (boost shared) pointer for this function class type
                                            SpudFunctionBucket_ptr;

}
#endif
