
#include "BoostTypes.h"
#include "SystemBucket.h"
#include "FunctionBucket.h"
#include "SolverBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Specific constructor
SystemBucket::SystemBucket()
{
  // Do nothing
}

// Specific constructor
SystemBucket::SystemBucket(Bucket* bucket) : bucket_(bucket)
{
  // Do nothing
}

// Default destructor
SystemBucket::~SystemBucket()
{
  empty_();
}

// Register a functionbucket as a field in the system
void SystemBucket::register_field(FunctionBucket_ptr field, std::string name)
{
  // First check if a field with this name already exists
  FunctionBucket_it f_it = fields_.find(name);
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

// Fetch a functionbucket as a field from the system
FunctionBucket_ptr SystemBucket::fetch_field(std::string name)
{
  // First check if a field with this name already exists
  FunctionBucket_it f_it = fields_.find(name);
  if (f_it == fields_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Field named \"%s\" does not exists in system.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

// Register a functionbucket as a coefficient in the system
void SystemBucket::register_coeff(FunctionBucket_ptr coeff, std::string name)
{
  // First check if a field with this name already exists
  FunctionBucket_it f_it = coeffs_.find(name);
  if (f_it != coeffs_.end())
  {
    // if it does, issue an error
    dolfin::error("Coefficient named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    coeffs_[name] = coeff;
  }
}

// Fetch a functionbucket as a coeff from the system
FunctionBucket_ptr SystemBucket::fetch_coeff(std::string name)
{
  // First check if a coeff with this name already exists
  FunctionBucket_it f_it = coeffs_.find(name);
  if (f_it == coeffs_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Coefficient named \"%s\" does not exists in system.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

// Register a ufl symbol-function name pair
void SystemBucket::register_uflname(std::string name, std::string uflsymbol)
{
  // First check if a name with this symbol already exists
  string_it s_it = uflnames_.find(uflsymbol);
  if (s_it != uflnames_.end())
  {
    // if it does, issue an error
    dolfin::error("Name with ufl symbol \"%s\" already exists in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then insert it into the maps
    uflnames_[uflsymbol] = name;
  }
}

// Return a pointer to a dolfin GenericFunction with the given uflsymbol
std::string SystemBucket::fetch_uflname(std::string uflsymbol)
{
  // First check if a function name with this symbol already exists
  string_it s_it = uflnames_.find(uflsymbol);
  if (s_it == uflnames_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Name with uflsymbol \"%s\" does not exist in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then return it
    return (*s_it).second;
  }
}

// Register a ufl symbol-function name pair
void SystemBucket::register_uflsymbol(GenericFunction_ptr function, std::string uflsymbol)
{
  // First check if a function with this symbol already exists
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);
  if (g_it != uflsymbols_.end())
  {
    // if it does, issue an error
    dolfin::error("GenericFunction with ufl symbol \"%s\" already exists in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then insert it into the maps
    uflsymbols_[uflsymbol] = function;
  }
}

// Create a ufl symbol-function name pair
void SystemBucket::create_uflsymbol(std::string uflsymbol)
{
  // First check if a function with this symbol already exists
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);
  if (g_it != uflsymbols_.end())
  {
    // if it does, issue an error
    dolfin::error("GenericFunction with ufl symbol \"%s\" already exists in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then insert it into the maps
    GenericFunction_ptr function;
    uflsymbols_[uflsymbol] = function;
  }
}

// Reset a ufl symbol-function name pair
void SystemBucket::reset_uflsymbol(GenericFunction_ptr function, std::string uflsymbol)
{
  // First check if a function with this symbol already exists
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);
  if (g_it == uflsymbols_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("GenericFunction with ufl symbol \"%s\" doesn't exist in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then insert it into the maps
    (*g_it).second = function;
  }
}

// Return a pointer to a dolfin GenericFunction with the given uflsymbol
GenericFunction_ptr SystemBucket::fetch_uflsymbol(std::string uflsymbol)
{
  // First check if a generic function with this symbol already exists
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);
  if (g_it == uflsymbols_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("GenericFunction with uflsymbol \"%s\" does not exist in system.", uflsymbol.c_str());
  }
  else
  {
    // if not then return it
    return (*g_it).second;
  }
}

// Register a coefficientspace in the system
void SystemBucket::register_coefficientspace(FunctionSpace_ptr coefficientspace, std::string name)
{
  // First check if a field with this name already exists
  FunctionSpace_it f_it = coefficientspaces_.find(name);
  if (f_it != coefficientspaces_.end())
  {
    // if it does, issue an error
    dolfin::error("FunctionSpace named \"%s\" already exists in system coefficientspaces.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    coefficientspaces_[name] = coefficientspace;
  }
}

// Check if the system contains a coefficientspace with the given name
bool SystemBucket::contains_coefficientspace(std::string name)
{
  // First check if a field with this name already exists
  FunctionSpace_it f_it = coefficientspaces_.find(name);
  return f_it != coefficientspaces_.end();
}

FunctionSpace_ptr SystemBucket::fetch_coefficientspace(std::string name)
{
  // First check if a functionspace with this name already exists
  FunctionSpace_it f_it = coefficientspaces_.find(name);
  if (f_it == coefficientspaces_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("FunctionSpace named \"%s\" doesn't exist in system coefficientspaces.", name.c_str());
  }
  else
  {
    // if it does return it
    return (*f_it).second;
  }
}

// Register a solver in the system
void SystemBucket::register_solver(SolverBucket_ptr solver, std::string name)
{
  // First check if a solver with this name already exists
  SolverBucket_it s_it = solvers_.find(name);
  if (s_it != solvers_.end())
  {
    // if it does, issue an error
    dolfin::error("SolverBucket named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    solvers_[name] = solver;
  }
}

// Return a string describing the contents of the system
std::string SystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << std::endl;
  indent++;
  s << uflsymbols_str(indent);
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  return s.str();
}

// Return a string describing the contents of uflsymbols_
std::string SystemBucket::uflsymbols_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( GenericFunction_const_it g_it = uflsymbols_.begin(); g_it != uflsymbols_.end(); g_it++ )
  {
    if ((*g_it).second)
    {
      s << indentation << "UFLSymbol " << (*g_it).first << " associated" << std::endl;
    }
    else
    {
      s << indentation << "UFLSymbol " << (*g_it).first << " not associated" << std::endl;
    }
  }
  return s.str();
}

// Return a string describing the contents of fields_
std::string SystemBucket::fields_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionBucket_const_it f_it = fields_.begin(); f_it != fields_.end(); f_it++ )
  {
    s << (*(*f_it).second).str(indent);
  }
  return s.str();
}

// Return a string describing the contents of fields_
std::string SystemBucket::coeffs_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionBucket_const_it f_it = coeffs_.begin(); f_it != coeffs_.end(); f_it++ )
  {
    s << (*(*f_it).second).str(indent);
  }
  return s.str();
}

// Return a string describing the contents of solvers_
std::string SystemBucket::solvers_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( SolverBucket_const_it s_it = solvers_.begin(); s_it != solvers_.end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
  }
  return s.str();
}

// Empty the system
void SystemBucket::empty_()
{
}

