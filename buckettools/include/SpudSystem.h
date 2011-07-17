
#ifndef __SPUD_SYSTEM_H
#define __SPUD_SYSTEM_H

#include "BoostTypes.h"
#include "System.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The SpudSystem class is a derived class of the system that populates the
  // data structures within a system using the spud options parser (and assumes
  // the structure of the buckettools schema)
  class SpudSystem : public System
  {
  // only accessible by this class
  private:

    // supplement the base class with an optionpath
    std::string optionpath_;

    // fill in some of the base info in the class
    void base_fill_();

    // fill in the system function information
    void systemfunction_fill_();

    // fill in the uflsymbols
    void uflsymbols_fill_();

    // fill in data about the fields
    void fields_fill_(uint &component,
                      std::map< uint, Expression_ptr > icexpressions);

    void apply_ic_(const uint &component, const std::map< uint, Expression_ptr > &icexpressions);

    // fill in data about the coefficients
    void coeffs_fill_();

  public:
    
    // Specific constructor
    SpudSystem(std::string optionpath, Bucket* bucket);

    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudSystem();

    // Fill the system assuming the buckettools schema
    void fill();

    // Return a string object describing the system
    std::string str() const
    { str(0); }

    // Return a string object describing the system
    std::string str(int indent) const;

    // Return the base optionpath for this system
    std::string optionpath() const
    { return optionpath_; }
    
  };
 
  // Define a boost shared pointer for a spud system
  typedef boost::shared_ptr< SpudSystem > SpudSystem_ptr;

}
#endif
