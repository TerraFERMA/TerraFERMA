
#include "DolfinBoostTypes.h"
#include "SpudFunctionBucket.h"
#include <dolfin.h>
#include <string>
#include <spud.h>
#include "SystemsWrapper.h"
#include "System.h"

using namespace buckettools;

// Specific constructor
SpudFunctionBucket::SpudFunctionBucket(std::string name, std::string optionpath, 
                                       GenericFunction_ptr function, System* system) : 
                                       optionpath_(optionpath), FunctionBucket(name, function, system)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudFunctionBucket::~SpudFunctionBucket()
{
  // Do nothing
}

// Fill the function using spud and assuming a buckettools schema structure
void SpudFunctionBucket::fill()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;

  buffer.str(""); buffer << optionpath() << "/type/output/include_in_diagnostics/functional";
  int nfuncs = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfuncs; i++)
  {
    buffer.str(""); buffer << optionpath() << "/type/output/include_in_diagnostics/functional[" << i << "]";
    functionals_fill_(buffer.str());
  }

}

void SpudFunctionBucket::functionals_fill_(const std::string &optionpath)
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
   
  std::string funcname;
  buffer.str(""); buffer << optionpath << "/name";
  Spud::get_option(buffer.str(), funcname);

  Form_ptr functional = ffc_fetch_functional((*system_).name(), name(), funcname, (*system_).mesh());

}

// Register a functional in the function
void SpudFunctionBucket::register_functional(Form_ptr functional, std::string name, std::string optionpath)
{
  // First check if a functional with this name already exists
  Form_it f_it = functionals_.find(name);
  if (f_it != functionals_.end())
  {
    // if it does, issue an error
    dolfin::error("Functional named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    functionals_[name]            = functional;
    functional_optionpaths_[name] = optionpath;
  }
}

// Return a string describing the contents of the spudfunction
std::string SpudFunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

