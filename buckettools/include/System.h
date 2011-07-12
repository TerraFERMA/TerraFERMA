
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  // Define iterator types for things accessed in the system maps (defined below)
  typedef std::map< std::string, FunctionSpace_ptr >::iterator        FunctionSpace_it;
  typedef std::map< std::string, FunctionSpace_ptr >::const_iterator  FunctionSpace_const_it;
  
  // The System class describes a functionspace and a set of solvers that act on the fields
  // contained in that (potentially mixed) functionspace.
  // This base class describes the basic data structures while derived classes may be defined
  // that allow it to be linked to an options system.
  class System
  {
  // only accessible by this class
  private:

    // the system name
    std::string name_;

    // Empty the data structures of the system
    void empty_();

  protected:
    
    // a pointer to the mesh which this system is on
    Mesh_ptr mesh_;

    // a pointer to the functionspace which this system describes
    FunctionSpace_ptr functionspace_;

    // the function (at several time-levels) on the above functionspace_
    Function_ptr function_, oldfunction_, iteratedfunction_;

    // a map from field names to the subfunctionspaces they are described on
    std::map< std::string, FunctionSpace_ptr > subfunctionspaces_;
    
  public:

    // No default constructor - always require a mesh pointer

    // Specific constructor with an uninitialised name
    System(Mesh_ptr mesh)
    { System("uninitialised_name", mesh); }

    // Specific constructor
    System(std::string name, Mesh_ptr mesh);
    
    // Default destructor
    ~System();

    // Return the system name
    std::string name()
    { return name_; }

    // Register a subfunctionspace in the system
    void register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name);

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< System > System_ptr;

}
#endif
