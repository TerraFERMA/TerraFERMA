
#include "DolfinBoostTypes.h"
#include "System.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Specific constructor
System::System(std::string name, Mesh_ptr mesh) : name_(name), mesh_(mesh)
{
  // Do nothing
}

// Default destructor
System::~System()
{
  empty_();
}

// Register a dolfin subfunctionspace in the system
void System::register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name)
{
  // First check if a subfunctionspace with this name already exists
  FunctionSpace_it f_it = subfunctionspaces_.find(name);
  if (f_it != subfunctionspaces_.end())
  {
    // if it does, issue an error
    dolfin::error("Subfunctionspace named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    subfunctionspaces_[name] = subfunctionspace;
  }
}

// Return whether a subfunctionspace with the given name is in the system
bool System::contains_subfunctionspace(std::string name)
{
  // First check if a subfunctionspace with this name already exists
  FunctionSpace_it f_it = subfunctionspaces_.find(name);
  return f_it != subfunctionspaces_.end();
}

// Return a pointer to a dolfin subfunctionspace with the given name
FunctionSpace_ptr System::fetch_subfunctionspace(std::string name)
{
  // First check if a subfunctionspace with this name already exists
  FunctionSpace_it f_it = subfunctionspaces_.find(name);
  if (f_it == subfunctionspaces_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Subfunctionspace named \"%s\" does not exist in system.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

// Register a dolfin expression for a bc in the system
void System::register_bcexpression(Expression_ptr bcexpression, std::string name)
{
  // First check if a bc expression with this name already exists
  Expression_it e_it = bcexpressions_.find(name);
  if (e_it != bcexpressions_.end())
  {
    // if it does, issue an error
    dolfin::error("BCExpression named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    bcexpressions_[name] = bcexpression;
  }
}

// Register a dolfin expression for a bc in the system
void System::register_dirichletbc(DirichletBC_ptr bc, std::string name)
{
  // First check if a bc with this name already exists
  DirichletBC_it bc_it = dirichletbcs_.find(name);
  if (bc_it != dirichletbcs_.end())
  {
    // if it does, issue an error
    dolfin::error("DirichletBC named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    dirichletbcs_[name] = bc;
  }
}

// Empty the spudbucket
void System::empty_()
{
  subfunctionspaces_.clear();
}

