
#include "BoostTypes.h"
#include "SpudFunctionBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemsWrapper.h"
#include "SpudSystem.h"
#include "SpudBase.h"

using namespace buckettools;

// Specific constructor
SpudFunctionBucket::SpudFunctionBucket(std::string name, std::string uflsymbol, std::string optionpath, 
                                       GenericFunction_ptr function, System* system) : 
                                       optionpath_(optionpath), FunctionBucket(name, uflsymbol, function, function, function, system)
{
  // Do nothing
}

// Specific constructor
SpudFunctionBucket::SpudFunctionBucket(std::string name, std::string uflsymbol, std::string optionpath, 
                                       GenericFunction_ptr function, GenericFunction_ptr oldfunction, 
                                       GenericFunction_ptr iteratedfunction, System* system) : 
                                       optionpath_(optionpath), FunctionBucket(name, uflsymbol, function, oldfunction, iteratedfunction, system)
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
  Spud::OptionError serr;
   
  std::string functionalname;
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), functionalname); spud_err(buffer.str(), serr);
  std::cout << "function name " << name() << std::endl;
  std::cout << "functional name " << functionalname << std::endl;

  Form_ptr functional = ufc_fetch_functional((*system_).name(), name(), functionalname, (*system_).mesh());
  register_functional(functional, functionalname, optionpath);

  // Loop over the functions requested by the form and hook up (potentially currently null) pointers
  uint ncoeff = (*functional).num_coefficients();
  for (uint i = 0; i < ncoeff; i++)
  {
    std::string uflsymbol = (*functional).coefficient_name(i);
    GenericFunction_ptr function = (*system_).fetch_uflsymbol(uflsymbol);
    if (!function)
    {
      std::string functionname = (*system_).fetch_uflname(uflsymbol);
      std::cout << "  requesting coefficient " << uflsymbol << ", " << functionname << std::endl; 
      std::cout << "  function not initialised" << std::endl;
      if (Spud::have_option((*dynamic_cast<SpudSystem*>(system_)).optionpath()+"/coefficient::"+functionname+"/type::Function"))
      {
        std::cout << "  but it is a function" << std::endl;
      }
      else
      {
        std::cout << "  and it is NOT a function" << std::endl;
      }

    }
    else
    {
      std::cout << "  requesting coefficient " << uflsymbol << std::endl; 
      std::cout << "  function is initialised" << std::endl;
    }

    (*functional).set_coefficient(uflsymbol, function);
  }

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

// Describe the contents of the functional_optionpaths_ map
std::string SpudFunctionBucket::functionals_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = functional_optionpaths_.begin(); s_it != functional_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Functional " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

