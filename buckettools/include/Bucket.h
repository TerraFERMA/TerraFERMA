
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
  // A class that describes a collection of meshes, systems, ufl symbols and point diagnostics (detectors).
  // The base class contains data structures and mechanisms for accessing those structures.
  // Derived classes are designed for specific options systems and are capable of populating those data
  // structures and may be extended to incorporate the ability to run entire models.
  //*****************************************************************|************************************************************//
  class Bucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    Bucket();                                                        // default constructor

    Bucket(std::string name);                                        // specific constructor specifying the name
    
    ~Bucket();                                                       // default destructor

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    virtual void run();                                              // run the model described in the bucket
 
    void solve();                                                    // solve all the systems in the bucket

    void update();                                                   // update the functions in the systems in this bucket

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a string containing the name of the bucket
    { return name_; }
    
    const int dimension() const                                      // return the geometry dimension
    { return dimension_; }

    //***************************************************************|***********************************************************//
    // Mesh data access
    //***************************************************************|***********************************************************//

    void register_mesh(Mesh_ptr mesh, std::string name);             // register a mesh with the given name in the bucket

    Mesh_ptr fetch_mesh(const std::string name);                     // return a (boost shared) pointer to a mesh with the given name
    
    Mesh_it meshes_begin();                                          // return an iterator to the beginning of the meshes

    Mesh_const_it meshes_begin() const;                              // return a constant iterator to the beginning of the meshes

    Mesh_it meshes_end();                                            // return an iterator to the end of the meshes

    Mesh_const_it meshes_end() const;                                // return a constant iterator to the end of the meshes
 
    //***************************************************************|***********************************************************//
    // System data access
    //***************************************************************|***********************************************************//

    void register_system(SystemBucket_ptr system, std::string name); // register a system with the given name in the bucket

    SystemBucket_ptr fetch_system(const std::string name);           // return a (boost shared) pointer to a system with the given name
    
    const SystemBucket_ptr fetch_system(const std::string name) 
          const;                                                     // return a constant (boost shared) pointer to a system with the given name
    
    SystemBucket_it systems_begin();                                 // return an iterator to the beginning of the systems

    SystemBucket_const_it systems_begin() const;                     // return a constant iterator to the beginning of the systems

    SystemBucket_it systems_end();                                   // return an iterator to the end of the systems

    SystemBucket_const_it systems_end() const;                       // return a constant iterator to the end of the systems

    int_SystemBucket_it orderedsystems_begin();                      // return an iterator to the beginning of the systems (in user
                                                                     //                                            prescribed order)

    int_SystemBucket_const_it orderedsystems_begin() const;          // return a constant iterator to the beginning of the systems (in user
                                                                     //                                            prescribed order)

    int_SystemBucket_it orderedsystems_end();                        // return an iterator to the end of the systems (in user
                                                                     //                                            prescribed order)

    int_SystemBucket_const_it orderedsystems_end() const;            // return a constant iterator to the end of the systems (in user
                                                                     //                                            prescribed order)

    //***************************************************************|***********************************************************//
    // UFL symbol data access
    //***************************************************************|***********************************************************//

    void register_baseuflsymbol(std::string name,                    // register a function name with a ufl symbol
                                std::string uflsymbol);  

    const std::string fetch_baseuflsymbol(std::string uflsymbol)     // return a string containing the function name belonging to a ufl symbol
                                                         const;

    const bool contains_baseuflsymbol(std::string uflsymbol) const;  // return a boolean, true if the bucket contains a ufl symbol,
                                                                     //                                            false otherwise

    void register_uflsymbol(GenericFunction_ptr function,            // register a (boost shared) pointer to a function with a ufl
                                   std::string uflsymbol);           // symbol

    GenericFunction_ptr fetch_uflsymbol(std::string uflsymbol)       // fetch a (boost shared) pointer to a function associated with
                       const;                                        // a ufl symbol

    //***************************************************************|***********************************************************//
    // Coefficient functionspace data access
    //***************************************************************|***********************************************************//

    void register_coefficientspace(FunctionSpace_ptr                 // register a functionspace for a coefficient function with a 
                           coefficientspace, std::string uflsymbol); // given ufl symbol

    const bool contains_coefficientspace(std::string uflsymbol)      // returns a boolean, true if the bucket contains a
                                                         const;      // functionspace for the uflsymbol, false otherwise

    FunctionSpace_ptr fetch_coefficientspace(std::string uflsymbol)  // returns a (boost shared) pointer to a functionspace for
                                                         const;      // the coefficient with the given uflsymbol

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    void output();                                                   // output diagnostics for the bucket

    const std::string str() const;                                   // return a string describing the bucket contents
    
    virtual const std::string meshes_str() const                     // return a string describing the meshes in the bucket
    { return meshes_str(0); }

    virtual const std::string meshes_str(int indent) const;          // return an indented string describing the meshes in the bucket

    const std::string systems_str() const                            // return a string describing the systems in the bucket
    { return systems_str(0); }

    const std::string systems_str(int indent) const;                 // return an indented string describing the systems in the bucket

    virtual const std::string coefficientspaces_str() const          // return a string describing the coefficient functionspaces
    { return coefficientspaces_str(0); }                             // contained in the bucket

    virtual const std::string coefficientspaces_str(int indent)      // return an indented string describing the coefficient functionspaces
                                                           const;    // contained in the bucket

    virtual const std::string uflsymbols_str() const                 // return a string describing the ufl symbols contained in the
    { return uflsymbols_str(0); }                                    // bucket

    virtual const std::string uflsymbols_str(int indent) const;      // return an indented string describing the ufl symbols contained in the
                                                                     // bucket

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to members of this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the name of the bucket
    
    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, std::string > baseuflsymbols_;            // a map from derived ufl symbols to their base symbol
    
    std::map< std::string, GenericFunction_ptr > uflsymbols_;        // a map from derived ufl symbols to field and coefficient pointers

    std::map< std::string, FunctionSpace_ptr > coefficientspaces_;   // a map from the base ufl symbol to a coefficient
                                                                     // functionspace

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the maps contained in the bucket

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // availble to this class and its derived classes

    //***************************************************************|***********************************************************//
    // Base data (continued)
    //***************************************************************|***********************************************************//

    int dimension_;                                                  // geometry dimension
                                                                     // (assumed size of various objects that don't state
                                                                     //  their own size explicitly)

    //***************************************************************|***********************************************************//
    // Pointers data (continued)
    //***************************************************************|***********************************************************//

    std::map< std::string, Mesh_ptr > meshes_;                       // a map from mesh names to (boost shared) pointers to meshes

    std::map< std::string, SystemBucket_ptr > systems_;              // a map from system names to (boost shared) pointers to systems

    std::map< int, SystemBucket_ptr > orderedsystems_;               // an ordered (user defined) map from system names to (boost
                                                                     // shared) pointers to systems

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void uflsymbols_fill_();                                         // fill the ufl symbol data structures

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

  typedef boost::shared_ptr< Bucket > Bucket_ptr;                    // define a boost shared ptr type for the class

}
#endif
