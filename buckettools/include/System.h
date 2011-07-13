
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  // Define iterator types for things accessed in the system maps (defined below)
  typedef std::map< std::string, FunctionSpace_ptr >::iterator       FunctionSpace_it;
  typedef std::map< std::string, FunctionSpace_ptr >::const_iterator FunctionSpace_const_it;
  typedef std::map< std::string, Function_ptr >::iterator            Function_it;
  typedef std::map< std::string, Function_ptr >::const_iterator      Function_const_it;
  typedef std::map< std::string, Expression_ptr >::iterator          Expression_it;
  typedef std::map< std::string, Expression_ptr >::const_iterator    Expression_const_it;
  typedef std::map< std::string, DirichletBC_ptr >::iterator         DirichletBC_it;
  typedef std::map< std::string, DirichletBC_ptr >::const_iterator   DirichletBC_const_it;
  typedef std::map< uint, Expression_ptr >::iterator                 uint_Expression_it;
  typedef std::map< uint, Expression_ptr >::const_iterator           uint_Expression_const_it;
  
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
    
    // a map from field names to the fields (subfunctions)
    std::map< std::string, Function_ptr > fields_;
    
    // a map from field::bc names to the subfunctionspaces they are described on
    std::map< std::string, Expression_ptr > bcexpressions_;
    
    // a map from field::bc::id names to dirichlet boundary conditions
    std::map< std::string, DirichletBC_ptr > dirichletbcs_;
    
    // a map from component integer index to initial condition expression
    std::map< uint, Expression_ptr > icexpressions_;
    
  public:

    // No default constructor - always require a mesh pointer

    // Specific constructor with an uninitialised name
    System(Mesh_ptr mesh)
    { System("uninitialised_name", mesh); }

    // Specific constructor
    System(std::string name, Mesh_ptr mesh);
    
    // Default destructor
    ~System();

    // Register a subfunctionspace in the system
    void register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name);

    // Return whether a subfunctionspace with the given name is in the system
    bool contains_subfunctionspace(std::string name);

    // Return a pointer to a dolfin subfunctionspace with the given name
    FunctionSpace_ptr fetch_subfunctionspace(std::string name);

    // Register a field (subfunction) in the system
    void register_field(Function_ptr field, std::string name);

    // Register a bc expression in the system
    void register_bcexpression(Expression_ptr bcexpression, std::string name);

    // Register a bc expression in the system
    void register_dirichletbc(DirichletBC_ptr bc, std::string name);

    // Register an initial condition expression
    void register_icexpression(Expression_ptr ic, uint component);

    // Return a string describing the contents of the system
    std::string str() const;

    // Print a description of the fields contained in the system
    virtual std::string fields_str() const;

    // Print a description of the bcexpressions contained in the system
    virtual std::string bcexpressions_str() const;

    // Return the system name
    std::string name() const
    { return name_; }

  };

  // Define a boost shared pointer to the system class type
  typedef boost::shared_ptr< System > System_ptr;

}
#endif
