
#ifndef __BUCKET_H
#define __BUCKET_H

// #include "GenericDetectors.h"
// #include "PythonDetectors.h"
#include "BoostTypes.h"
#include "System.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The Bucket is a class that describes a collection of meshes, systems and point diagnostics (detectors).
  // The base class contains data structures and mechanisms for accessing those structures.
  // Derived classes are designed for specific options systems and are capable of populating those data
  // structures and may be extended to incorporate the ability to run entire models.
  class Bucket
  {
  // only accessible to memebers of this class
  private:

    // The name of the bucket
    std::string name_;
    
    // Destroy the data in the bucket
    void empty_();

  // accessible to derived classes
  protected:    

    // Geometry dimension (assumed size of various objects that don't state their own size explicitly)
    int dimension_;

    // A map from mesh name to mesh pointer
    // The mesh format is assumed to be DOLFIN
    std::map< std::string, Mesh_ptr >         meshes_;

    // A map from system name to system pointer
    // Systems describe function spaces and the solvers that operate on them
    std::map< std::string, System_ptr >       systems_;

//    // A map from detector set name to detectors
//    std::map< std::string, GenericDetectors_ptr >    detectors_;
    
  // accessible to everyone
  public:
    
    // Default constructor (the name really isn't all that necessary)
    Bucket();

    // Optional constructor specifying the name
    Bucket(std::string name);
    
    // Default destructor
    ~Bucket();

    // Print a description of the bucket contents
    std::string str() const;
    
    // Print a description of the meshes contained in the bucket
    virtual std::string meshes_str() const
    { meshes_str(0); }

    // Print a description of the meshes contained in the bucket
    virtual std::string meshes_str(int indent) const;

    // Print a description of the systems contained in the bucket
    std::string systems_str() const
    { systems_str(0); }

    // Print a description of the systems contained in the bucket
    std::string systems_str(int indent) const;

    // Print the name of the bucket
    std::string name() const
    { return name_; }
    
    // Tools for manipulating the base bucket class data structures:
    // Register a mesh in the bucket (i.e. put it into the meshes_ map)
    void register_mesh(Mesh_ptr mesh, std::string name);

    // Return a mesh pointer from the meshes_ map, given its name
    Mesh_ptr fetch_mesh(const std::string name);
    
    // Return the iterator to the beginning of the meshes_ map
    Mesh_it meshes_begin();

    // Return the const iterator to the beginning of the meshes_ map
    Mesh_const_it meshes_begin() const;

    // Return the iterator to the end of the meshes_ map
    Mesh_it meshes_end();

    // Return the const iterator to the end of the meshes_ map
    Mesh_const_it meshes_end() const;
 
    // Register a system in the bucket (i.e. put it into the systems_ map)
    void register_system(System_ptr system, std::string name);

    // Return a system pointer from the systems_ map, given its name
    System_ptr fetch_system(const std::string name);
    
    // Return the iterator to the beginning of the systems_ map
    System_it systems_begin();

    // Return the const iterator to the beginning of the systems_ map
    System_const_it systems_begin() const;

    // Return the iterator to the end of the systems_ map
    System_it systems_end();

    // Return the const iterator to the end of the systems_ map
    System_const_it systems_end() const;

    // Return the geometry dimension
    int dimension() const
    { return dimension_; }
 
//    void register_system(System_ptr std::string name)
//    { register_system(name, "uninitialised_path"); }
//
//    void register_system(std::string name, std::string option_path);
//
//    void register_detector(GenericDetectors_ptr detector, std::string name)
//    { register_detector(detector, name, "uninitialised_path"); }
//    
//    void register_detector(GenericDetectors_ptr detector, std::string name, std::string option_path);
//    
//    void register_functionspace(FunctionSpace_ptr functionspace, std::string systemname, std::string name);
//    
//    void register_dirichletbc(DirichletBC_ptr dirichletbc, std::string systemname, std::string name);
//    
//    void register_function(Function_ptr function, std::string systemname, std::string name)
//    { register_function(function, systemname, name, "uninitialised_path"); }
//    
//    void register_function(Function_ptr function, std::string systemname, std::string name, std::string option_path);
//    
//    FunctionSpace_ptr fetch_functionspace(const std::string systemname, const std::string name);
//    
//    Function_ptr fetch_function(const std::string systemname, const std::string name);
//    
//    GenericDetectors_ptr fetch_detector(const std::string name);
//    
//
//    std::map< std::string, GenericDetectors_ptr >::iterator detectors_begin();
//    
//    std::map< std::string, GenericDetectors_ptr >::const_iterator detectors_begin() const;
//    
//    std::map< std::string, GenericDetectors_ptr >::iterator detectors_end();
//    
//    std::map< std::string, GenericDetectors_ptr >::const_iterator detectors_end() const;
//    
//    std::map< std::string, DirichletBC_ptr >::iterator dirichletbcs_begin();
//    
//    std::map< std::string, DirichletBC_ptr >::const_iterator dirichletbcs_begin() const;
//    
//    std::map< std::string, DirichletBC_ptr >::iterator dirichletbcs_end();
//    
//    std::map< std::string, DirichletBC_ptr >::const_iterator dirichletbcs_end() const;
//    
//    std::map< std::string, Function_ptr >::iterator functions_begin();
//    
//    std::map< std::string, Function_ptr >::const_iterator functions_begin() const;
//    
//    std::map< std::string, Function_ptr >::iterator functions_end();
//    
//    std::map< std::string, Function_ptr >::const_iterator functions_end() const;
//    
//    void register_bcexp(GenericFunction_ptr bcexp, std::string systemname);
    
  };

  // Define a boost shared ptr type for the class
  typedef boost::shared_ptr< Bucket > Bucket_ptr;

}
#endif
