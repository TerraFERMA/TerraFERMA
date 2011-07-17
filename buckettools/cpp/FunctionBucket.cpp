
#include "BoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
FunctionBucket::FunctionBucket()
{
  // Do nothing
}

// Specific constructor
FunctionBucket::FunctionBucket(System* system) : system_(system)
{
  // Do nothing
}

// Default destructor
FunctionBucket::~FunctionBucket()
{
  empty_();
}

// Register a dolfin subfunctionspace in the function
void FunctionBucket::register_subfunctionspace(FunctionSpace_ptr subfunctionspace, int index)
{
  // First check if a subfunctionspace with this name already exists
  int_FunctionSpace_it f_it = subfunctionspaces_.find(index);
  if (f_it != subfunctionspaces_.end())
  {
    // if it does, issue an error
    dolfin::error("Subfunctionspace with index \"%d\" already exists in system.", index);
  }
  else
  {
    // if not then insert it into the maps
    subfunctionspaces_[index] = subfunctionspace;
  }
}

// Return whether a subfunctionspace with the given component index is in the system
bool FunctionBucket::contains_subfunctionspace(int index)
{
  // First check if a subfunctionspace with this name already exists
  int_FunctionSpace_it f_it = subfunctionspaces_.find(index);
  return f_it != subfunctionspaces_.end();
}

// Return a pointer to a dolfin subfunctionspace with the given component index
FunctionSpace_ptr FunctionBucket::fetch_subfunctionspace(int index)
{
  // First check if a subfunctionspace with this name already exists
  int_FunctionSpace_it f_it = subfunctionspaces_.find(index);
  if (f_it == subfunctionspaces_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Subfunctionspace with index \"%d\" does not exist in system.", index);
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

// Register a functional in the function
void FunctionBucket::register_functional(Form_ptr functional, std::string name)
{
  // First check if a field with this name already exists
  Form_it f_it = functionals_.find(name);
  if (f_it != functionals_.end())
  {
    // if it does, issue an error
    dolfin::error("Functional named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    functionals_[name] = functional;
  }
}

// Return a pointer to a functional with the given name
Form_ptr FunctionBucket::fetch_functional(std::string name)
{
  // First check if a functional with this name already exists
  Form_it f_it = functionals_.find(name);
  if (f_it == functionals_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Functional named \"%s\" does not exist in function.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

// Register a dolfin expression for a bc in the functionbucket
void FunctionBucket::register_bcexpression(Expression_ptr bcexpression, std::string name)
{
  // First check if a bc expression with this name already exists
  Expression_it e_it = bcexpressions_.find(name);
  if (e_it != bcexpressions_.end())
  {
    // if it does, issue an error
    dolfin::error("BCExpression named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    bcexpressions_[name] = bcexpression;
  }
}

// Register a dolfin expression for a bc in the functionbucket
void FunctionBucket::register_dirichletbc(DirichletBC_ptr bc, std::string name)
{
  // First check if a bc with this name already exists
  DirichletBC_it bc_it = dirichletbcs_.find(name);
  if (bc_it != dirichletbcs_.end())
  {
    // if it does, issue an error
    dolfin::error("DirichletBC named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    dirichletbcs_[name] = bc;
  }
}

// Return a string describing the contents of the function
std::string FunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

// Return a string describing the contents of functionals_
std::string FunctionBucket::functionals_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = functionals_.begin(); f_it != functionals_.end(); f_it++ )
  {
    s << indentation << "Functional " << (*f_it).first  << std::endl;
  }
  return s.str();
}

// Empty the function
void FunctionBucket::empty_()
{
  functionals_.clear();
}

