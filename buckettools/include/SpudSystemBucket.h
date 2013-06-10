
#ifndef __SPUD_SYSTEM_H
#define __SPUD_SYSTEM_H

#include "BoostTypes.h"
#include "SystemBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudSystemBucket class:
  //
  // The SpudSystemBucket class is a derived class of the system that populates the
  // data structures within a system using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudSystemBucket : public SystemBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudSystemBucket(const std::string &optionpath, Bucket* bucket); // default constructor

    ~SpudSystemBucket();                                             // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill the system assuming a buckettools spud schema

    void allocate_coeff_function();                                  // allocate the coefficient functions

    void initialize_forms();                                         // initialize the forms and functionals by attaching the coeffs

    void initialize_fields_and_coefficients();                       // initialize the expressions that fields and coefficients use

    void initialize_solvers();                                       // initialize the matrices related to forms as well as the
                                                                     // petsc objects

    void copy_diagnostics(SystemBucket_ptr &system, 
                            Bucket_ptr &bucket) const;               // copy the data necessary for the diagnostics data file(s)

    std::vector< GenericFunction_ptr > collect_vis_functions() const;// output the diagnostics on this system

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a string containing the optionpath for the system
    { return optionpath_; }
    
    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str() const                                    // return a string describing the contents of the system
    { return str(0); }

    const std::string str(int indent) const;                         // return an indented string describing the contents of the
                                                                     // system

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible in this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the system optionpath

    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill in the base information about the system

    void fill_systemfunction_();                                     // fill in the system function information

    void fill_fields_();                                             // fill in the data about the system fields (subfunctions)

    void fill_bcs_();                                                // fill in the data about the system bcs

    void fill_points_();                                             // fill in the data about the system points

    void fill_coeffs_();                                             // fill in the coefficient information

    void fill_solvers_();                                            // fill in the solver bucket information

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    void checkpoint_options_();                                      // checkpoint the options system for the spudsystembucket

  };
 
  typedef boost::shared_ptr< SpudSystemBucket > SpudSystemBucket_ptr;// define a (boost shared) pointer type for the spud system

}
#endif
