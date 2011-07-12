
#ifndef __SPUD_BUCKET_H
#define __SPUD_BUCKET_H

#include "Bucket.h"
#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The SpudBucket derived class uses the data structures of the Bucket class,
  // which describes a numerical model, and populates them (plus some supplementary
  // material) using the spud options parser.
  class SpudBucket : public Bucket
  {
  // only accessible to this class
  private:
    
    // the spud xml base option path for this bucket (probably an empty string)
    std::string optionpath_;

    // supplementary to the bucket base data store the optionpaths for the meshes_
    std::map< std::string, std::string > mesh_optionpaths_;

    // supplementary to the bucket base data store the optionpaths for the systems_
    // (superfluous? as each spudsystem derived class can also store this info, unlike meshes)
    std::map< std::string, std::string > system_optionpaths_;

    // populate the meshes_ data structure
    void meshes_fill_(const std::string &optionpath);
    
    // populate the systems_ data structure
    void systems_fill_(const std::string &optionpath);

//    // populate the detectors_ data structure
//    void detectors_fill_();
 
    // return the optionpath for this bucket
    std::string optionpath()
    { return optionpath_; }

  public:
    
    // Default constructor
    SpudBucket()
    { SpudBucket("uninitialised_name", ""); }
    
    // Specific constructor (assume string is a name and optionpath is empty)
    SpudBucket(std::string name)
    { SpudBucket(name, ""); }
    
    // Specific constructor
    SpudBucket(std::string name, std::string option_path);

    // Destructor (virtual so that it calls the base destructor as well)
    virtual ~SpudBucket();
    
    // Populate the bucket assuming the buckettools spud schema
    void fill();

    // Tools for manipulating the class data structures:
    // Register a mesh in the bucket with an optionpath in the derived class
    void register_mesh(Mesh_ptr mesh, std::string name, std::string optionpath);

    // Return the optionpath to a named mesh
    std::string fetch_mesh_optionpath(const std::string name);

    // Describe the meshes in the spudbucket (including optionpaths)
    std::string meshes_str() const;

    // Return the iterator to the beginning of the mesh_optionpaths_ map
    string_it mesh_optionpaths_begin();

    // Return the const iterator to the beginning of the mesh_optionpaths_ map
    string_const_it mesh_optionpaths_begin() const;

    // Return the iterator to the end of the mesh_optionpaths_ map
    string_it mesh_optionpaths_end();

    // Return the const iterator to the end of the mesh_optionpaths_ map
    string_const_it mesh_optionpaths_end() const;
 
  };
}
#endif
