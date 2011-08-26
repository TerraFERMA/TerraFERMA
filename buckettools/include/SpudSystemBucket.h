
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

    SpudSystemBucket(std::string optionpath, Bucket* bucket);        // default constructor

    ~SpudSystemBucket();                                     // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill the system assuming a buckettools spud schema

    void funccoeffs_fill();                                          // fill in the coefficient functions

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

    void base_fill_();                                               // fill in the base information about the system

    void systemfunction_fill_();                                     // fill in the system function information

    void fields_fill_();                                             // fill in the data about the system fields (subfunctions)

    void apply_ic_(const uint &component, const std::map< uint,      // apply the initial conditions to the system function
                                  Expression_ptr > &icexpressions);

    void apply_bc_();                                                // apply the bcs to the system function

    void expcoeffs_fill_();                                          // fill in the coefficient expression information

    void solvers_fill_();                                            // fill in the solver bucket information

  };
 
  typedef boost::shared_ptr< SpudSystemBucket > SpudSystemBucket_ptr;// define a (boost shared) pointer type for the spud system

}
#endif
