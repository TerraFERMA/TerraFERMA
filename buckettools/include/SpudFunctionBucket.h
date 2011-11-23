
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

    SpudFunctionBucket(const std::string &optionpath, 
                                               SystemBucket* system);// specific constructor
    
    ~SpudFunctionBucket();                                           // default destructor

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a string containing the optionpath of this function
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill_field(const uint &index);                              // fill this function assuming it's a field

    void fill_coeff(const uint &index);                              // fill this function assuming it's a coefficient

    void allocate_coeff_function();                                  // allocate this function assuming it's a coefficient function

    void initialize_field();                                         // initialize the expressions associated with a field

    void initialize_coeff_expression();                              // initialize the expressions associated with a field

    void initialize_coeff_function();                                // initialize the expressions associated with a field

    void copy_diagnostics(FunctionBucket_ptr &function, 
                                  SystemBucket_ptr &system) const;   // copy the data necessary for the diagnostics data file(s)

    //***************************************************************|***********************************************************//
    // Functional data access
    //***************************************************************|***********************************************************//

    void register_functional(Form_ptr functional, 
                                    const std::string &name,         // register a functional with the given name and optionpath
                                    const std::string &optionpath);

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const bool include_in_visualization() const;                     // return a boolean indicating if this function is to 
                                                                     // be included in diagnostic output
    
    const bool include_in_statistics() const;                        // return a boolean indicating if this function is to 
                                                                     // be included in diagnostic output
    
    const bool include_in_steadystate() const;                       // return a boolean indicating if this function is to 
                                                                     // be included in steadystate output
    
    const bool include_in_detectors() const;                         // return a boolean indicating if this function is to 
                                                                     // be included in steadystate output
    
    const std::string str() const                                    // return a string describing the contents of this function
    { return str(0); }

    const std::string str(int indent) const;                         // return an indented string describing the contents of this function

    const std::string functionals_str() const                        // return a string describing the functionals of this function
    { return str(0); }

    const std::string functionals_str(const int &indent) const;      // return an indented string describing the functionals of this function

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

    void fill_base_(const uint &index);                              // fill the base data of this function
 
    void allocate_field_();                                          // allocate this function assuming its a field

    void allocate_coeff_expression_();                               // allocate this function assuming its a coefficient expression

    void fill_bc_(const std::string &optionpath);                    // fill in the bc for this function

    void fill_bc_component_(const std::string &optionpath,           // fill in the bc data for a component of this function
                            const std::string &bcname,
                            const std::vector<int> &bcids);

    void initialize_bc_(const std::string &optionpath);              // fill in the bc for this function

    void initialize_bc_component_(const std::string &optionpath,     // fill in the bc data for a component of this function
                                  const std::string &bcname);

    void fill_constantfunctional_();                                 // fill in the data for constant expressions defined by functionals

    void fill_functionals_();                                        // fill in the data for the functionals of this function

    Expression_ptr allocate_expression_over_regions_(
                                      const std::string &optionpath,
                                      const double_ptr time);        // allocate an expression over region ids based on an optionpath

    Expression_ptr allocate_expression_over_regions_(
                                      const std::string &optionpath,
                                      const double_ptr time,
                                      bool *time_dependent);         // allocate an expression over region ids based on an optionpath

    Expression_ptr allocate_expression_(
                                      const std::string &optionpath,
                                      const std::string &expressionname,
                                      const double_ptr time);        // allocate an expression based on an optionpath

    Expression_ptr allocate_expression_(
                                      const std::string &optionpath,
                                      const std::string &expressionname,
                                      const double_ptr time,
                                      bool *time_dependent);         // allocate an expression based on an optionpath

    void initialize_expression_over_regions_(
                                      Expression_ptr expression,
                                      const std::string &optionpath);// allocate an expression over region ids based on an optionpath

    void initialize_expression_(Expression_ptr expression,
                                const std::string &optionpath,
                                const std::string &expressionname);  // allocate an expression based on an optionpath

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                           // empty the data structures of this function bucket

  };
 
  typedef boost::shared_ptr< SpudFunctionBucket >                    // define a (boost shared) pointer for this function class type
                                            SpudFunctionBucket_ptr;

}
#endif
