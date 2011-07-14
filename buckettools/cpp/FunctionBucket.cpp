
#include "DolfinBoostTypes.h"
#include "FunctionBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Specific constructor
FunctionBucket::FunctionBucket(std::string name, GenericFunction_ptr function, System* system) : name_(name), function_(function), system_(system)
{
  // Do nothing
}

// Default destructor
FunctionBucket::~FunctionBucket()
{
  empty_();
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

