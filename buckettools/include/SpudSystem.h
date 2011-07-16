
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

    // supplementary to the system base data store the optionpaths for the fields_
    std::map< std::string, std::string > field_optionpaths_;

    // fill in the uflsymbols
    void uflsymbols_fill_();

    // fill in data about the fields
    void fields_fill_(const std::string &optionpath, 
                      const uint &field_i, 
                      const uint &nfields,
                      uint &component);

    // fill in data about the coefficients
    void coeffs_fill_(const std::string &optionpath);

    // fill in data about the bcs
    void bc_fill_(const std::string &optionpath,
                  const std::string &fieldname,
                  const int &size,
                  const std::vector<int> &shape,
                  const FunctionSpace_ptr &subfunctionspace, 
                  const MeshFunction_uint_ptr &edgeidmeshfunction);

    // fill in data about a bc component
    void bc_component_fill_(const std::string &optionpath,
                            const std::string &fieldname,
                            const std::string &bcname,
                            const int &size,
                            const std::vector<int> &shape,
                            const std::vector<int> &bcids,
                            const FunctionSpace_ptr &subfunctionspace,
                            const MeshFunction_uint_ptr &edgeidmeshfunction);

    // fill in data about an individual field's ic
    void ic_fill_(const std::string &optionpath,
                  const int &size,
                  const std::vector<int> &shape,
                  uint &component);

  public:
    
    // FIXME: don't know why this constructor doesn't override base class constructor
    // Specific constructor with no name or path
    //SpudSystem(Mesh_ptr mesh)
    //{ SpudSystem("uninitialised_name", "uninitialised_path", mesh); }

    // FIXME: don't know why this constructor doesn't override base class constructor
    // Specific constructor assuming no path
    //SpudSystem(std::string name, Mesh_ptr mesh)
    //{ SpudSystem(name, "uninitialised_path", mesh); }

    // Specific constructor
    SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh, Bucket* bucket);
    
    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudSystem();

    // Fill the system assuming the buckettools schema
    void fill();

    // Register a dolfin function as a field in the system (with a spud optionpath)
    void register_field(FunctionBucket_ptr field, std::string name, std::string optionpath);

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
