
#ifndef __SPUD_SYSTEM_H
#define __SPUD_SYSTEM_H

#include "BoostTypes.h"
#include "SystemBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The SpudSystemBucket class is a derived class of the system that populates the
  // data structures within a system using the spud options parser (and assumes
  // the structure of the buckettools schema)
  class SpudSystemBucket : public SystemBucket
  {
  // only accessible by this class
  private:

    // supplement the base class with an optionpath
    std::string optionpath_;

    // fill in some of the base info in the class
    void base_fill_();

    // fill in the system function information
    void systemfunction_fill_();

    // fill in data about the fields
    void fields_fill_();

    void apply_ic_(const uint &component, const std::map< uint, Expression_ptr > &icexpressions);

    void apply_bc_();

    // fill in data about the coefficients
    void expcoeffs_fill_();

    // fill in data about the solvers
    void solvers_fill_();

  public:
    
    // Specific constructor
    SpudSystemBucket(std::string optionpath, Bucket* bucket);

    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudSystemBucket();

    // Fill the system assuming the buckettools schema
    void fill();

    // Fill the aliased components of the system assuming the buckettools schema
    void attach_and_initialize();

    // fill in data about the coefficients
    void funccoeffs_fill();

    // Return a string object describing the system
    const std::string str() const
    { return str(0); }

    // Return a string object describing the system
    const std::string str(int indent) const;

    // Return the base optionpath for this system
    const std::string optionpath() const
    { return optionpath_; }
    
  };
 
  // Define a boost shared pointer for a spud system
  typedef boost::shared_ptr< SpudSystemBucket > SpudSystemBucket_ptr;

}
#endif
