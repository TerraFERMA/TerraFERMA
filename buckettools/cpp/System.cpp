
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

// Register a dolfin function as a field in the system
void System::register_field(Function_ptr field, std::string name)
{
  // First check if a field with this name already exists
  Function_it f_it = fields_.find(name);
  if (f_it != fields_.end())
  {
    // if it does, issue an error
    dolfin::error("Field named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    fields_[name] = field;
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

// Register a dolfin expression for an ic in the system
void System::register_icexpression(Expression_ptr ic, uint component)
{
  // First check if an ic with this component index already exists
  uint_Expression_it ic_it = icexpressions_.find(component);
  if (ic_it != icexpressions_.end())
  {
    // if it does, issue an error
    dolfin::error("ICExpression with component index \"%d\" already exists in system.", component);
  }
  else
  {
    // if not then insert it into the maps
    icexpressions_[component] = ic;
  }
}

// Empty the spudbucket
void System::empty_()
{
  subfunctionspaces_.clear();
}

