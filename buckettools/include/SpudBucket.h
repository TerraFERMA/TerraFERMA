
#ifndef __SPUD_BUCKET_H
#define __SPUD_BUCKET_H

#include "Bucket.h"
#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudBucket class:
  //
  // The SpudBucket derived class uses the data structures of the Bucket class,
  // which describes a numerical model, and populates them (plus some supplementary
  // material) using the spud options parser.
  //*****************************************************************|************************************************************//
  class SpudBucket : public Bucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudBucket();                                                    // Default constructor
    
    SpudBucket(std::string name);                                    // Optional constructor (specifying name)
    
    SpudBucket(std::string name, std::string option_path);           // Optional constructor (specifying name and optionpath)

    ~SpudBucket();                                           // Default destructor
    
    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    void run();                                                      // run the model described in the bucket
 
    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill the bucket data structures assuming the buckettools
                                                                     // spud schema

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath()                                   // return the optionpath for this bucket
    { return optionpath_; }                                          // (normally an empty string)

    //***************************************************************|***********************************************************//
    // Mesh data access
    //***************************************************************|***********************************************************//

    void register_mesh(Mesh_ptr mesh, std::string name,              // register a mesh with a given name (and an optionpath)
                                        std::string optionpath); 

    std::string fetch_mesh_optionpath(const std::string name);       // return the optionpath associated with the named mesh

    string_it mesh_optionpaths_begin();                              // return an iterator to the beginning of the mesh optionpaths

    string_const_it mesh_optionpaths_begin() const;                  // return a constant iterator to the beginning of the mesh
                                                                     // optionpaths

    string_it mesh_optionpaths_end();                                // return an iterator to the end of the mesh optionpaths

    string_const_it mesh_optionpaths_end() const;                    // return a constant iterator to the end of the mesh
                                                                     // optionpaths
 
    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string meshes_str() const                             // return a string describing the meshes
    { return meshes_str(0); }

    const std::string meshes_str(int indent) const;                  // return an indented string describing the meshes

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to this class
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the option path associated with the bucket (normally an
                                                                     // empty string)

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, std::string > mesh_optionpaths_;          // a map from mesh names to spud mesh optionpaths

    //***************************************************************|***********************************************************//
    // Functions used to run the model (continued)
    //***************************************************************|***********************************************************//

    void timestep_run_();                                            // run a dynamics simulation

    void steady_run_();                                              // run a steady state simulation

    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void meshes_fill_(const std::string &optionpath);                // fill in the mesh data structures
    
    void systems_fill_(const std::string &optionpath);               // fill in information about the systems

    void baseuflsymbols_fill_(const std::string &optionpath);        // fill the ufl symbol maps

//    void detectors_fill_();                                          // fill the detectors
 
    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void derived_empty_();

  };

  typedef boost::shared_ptr< SpudBucket > SpudBucket_ptr;            // define a boost shared pointer type for this class

}
#endif
