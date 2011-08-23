
#ifndef __BUCKET_H
#define __BUCKET_H

// #include "GenericDetectors.h"
// #include "PythonDetectors.h"
#include "BoostTypes.h"
#include "SystemBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // Bucket class:
  //
  // A class that describes a collection of meshes, systems and point diagnostics (detectors).
  // The base class contains data structures and mechanisms for accessing those structures.
  // Derived classes are designed for specific options systems and are capable of populating those data
  // structures and may be extended to incorporate the ability to run entire models.
  //*****************************************************************|************************************************************//
  class Bucket
  {
  public:                                                            // accessible to everyone
    
    Bucket();                                                        // default constructor

    Bucket(std::string name);                                        // optional constructor specifying the name
    
    ~Bucket();                                                       // default destructor

    const std::string str() const                                    // return a string describing the bucket contents
    
    virtual const std::string meshes_str() const                     // return a string describing the meshes in the bucket
    { return meshes_str(0); }

    virtual const std::string meshes_str(int indent) const;          // return an indented string describing the meshes in the bucket

    const std::string systems_str() const                            // return a string describing the systems in the bucket
    { return systems_str(0); }

    const std::string systems_str(int indent) const;                 // return an indented string describing the systems in the bucket

    const std::string name() const                                   // return a string containing the name of the bucket
    { return name_; }
    
    void register_mesh(Mesh_ptr mesh, std::string name);             // register a mesh with the given name in the bucket

    Mesh_ptr fetch_mesh(const std::string name);                     // return a (boost shared) pointer to a mesh with the given name
    
    Mesh_it meshes_begin();                                          // return an iterator to the beginning of the meshes

    Mesh_const_it meshes_begin() const;                              // return a constant iterator to the beginning of the meshes

    Mesh_it meshes_end();                                            // return an iterator to the end of the meshes

    Mesh_const_it meshes_end() const;                                // return a constant iterator to the end of the meshes
 
    void register_system(SystemBucket_ptr system, std::string name); // register a system with the given name in the bucket

    SystemBucket_ptr fetch_system(const std::string name);           // return a (boost shared) pointer to a system with the given name
    
    const SystemBucket_ptr fetch_system(const std::string name) 
          const;                                                     // return a constant (boost shared) pointer to a system with the given name
    
    SystemBucket_it systems_begin();

    // Return the const iterator to the beginning of the systems_ map
    SystemBucket_const_it systems_begin() const;

    // Return the iterator to the end of the systems_ map
    SystemBucket_it systems_end();

    // Return the const iterator to the end of the systems_ map
    SystemBucket_const_it systems_end() const;

    // Return the iterator to the beginning of the systems_ map
    int_SystemBucket_it orderedsystems_begin();

    // Return the const iterator to the beginning of the systems_ map
    int_SystemBucket_const_it orderedsystems_begin() const;

    // Return the iterator to the end of the systems_ map
    int_SystemBucket_it orderedsystems_end();

    // Return the const iterator to the end of the systems_ map
    int_SystemBucket_const_it orderedsystems_end() const;

    // Register a ufl name
    void register_uflname(std::string name, std::string uflsymbol);

    // Return the name of a function with uflsymbol
    const std::string fetch_uflname(std::string uflsymbol) const;

    // Return the name of a function with uflsymbol
    const bool contains_uflname(std::string uflsymbol) const;

    // Register a ufl system and function pointer in the system
    void register_uflsymbol(GenericFunction_ptr function, std::string uflsymbol);

    // Return a pointer to a dolfin GenericFunction with the given uflsymbol
    GenericFunction_ptr fetch_uflsymbol(std::string uflsymbol) const;

    void register_coefficientspace(FunctionSpace_ptr coefficientspace, std::string name);

    const bool contains_coefficientspace(std::string name) const;

    FunctionSpace_ptr fetch_coefficientspace(std::string name) const;

    void solve();

    void update();

    // Return the geometry dimension
    const int dimension() const
    { return dimension_; }

    void output();

    virtual void run();
 
    // Print a description of the coefficient functionspaces contained in the system
    virtual const std::string coefficientspaces_str() const
    { return coefficientspaces_str(0); }

    // Print a description of the coefficient functionspaces contained in the system
    virtual const std::string coefficientspaces_str(int indent) const;

    // Print a description of the fields contained in the system
    virtual const std::string uflsymbols_str() const
    { return uflsymbols_str(0); }

    // Print a description of the fields contained in the system
    virtual const std::string uflsymbols_str(int indent) const;

  // only accessible to memebers of this class
  private:

    std::string name_;                                                             // The name of the bucket
    
    // a map from ufl symbols to field and coefficient names
    std::map< std::string, std::string > uflnames_;
    
    // a map from ufl symbols to pointers to functions
    std::map< std::string, GenericFunction_ptr > uflsymbols_;

    std::map< std::string, FunctionSpace_ptr > coefficientspaces_;

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
    // SystemBuckets describe function spaces and the solvers that operate on them
    std::map< std::string, SystemBucket_ptr >       systems_;
    std::map< int, SystemBucket_ptr >               orderedsystems_;

    // fill in the uflsymbols
    void uflsymbols_fill_();

//    // A map from detector set name to detectors
//    std::map< std::string, GenericDetectors_ptr >    detectors_;
    
//    void register_system(SystemBucket_ptr std::string name)
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
