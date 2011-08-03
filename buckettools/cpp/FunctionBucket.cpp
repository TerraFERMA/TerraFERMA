
#include "BoostTypes.h"
#include "FunctionBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
FunctionBucket::FunctionBucket()
{
  // Do nothing
}

// Specific constructor
FunctionBucket::FunctionBucket(SystemBucket* system) : system_(system)
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

Form_it FunctionBucket::functionals_begin()
{
  return functionals_.begin();
}

Form_const_it FunctionBucket::functionals_begin() const
{
  return functionals_.begin();
}

Form_it FunctionBucket::functionals_end()
{
  return functionals_.end();
}

Form_const_it FunctionBucket::functionals_end() const
{
  return functionals_.end();
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
void FunctionBucket::register_bc(BoundaryCondition_ptr bc, std::string name)
{
  // First check if a bc with this name already exists
  BoundaryCondition_it bc_it = bcs_.find(name);
  if (bc_it != bcs_.end())
  {
    // if it does, issue an error
    dolfin::error("BoundaryCondition named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    bcs_[name] = bc;
  }
}

BoundaryCondition_it FunctionBucket::bcs_begin()
{
  return bcs_.begin();
}

BoundaryCondition_const_it FunctionBucket::bcs_begin() const
{
  return bcs_.begin();
}

BoundaryCondition_it FunctionBucket::bcs_end()
{
  return bcs_.end();
}

BoundaryCondition_const_it FunctionBucket::bcs_end() const
{
  return bcs_.end();
}

void FunctionBucket::attach_functional_coeffs()
{
  for (Form_it f_it = functionals_begin(); f_it != functionals_end(); f_it++)
  {
    dolfin::info("  Attaching coeffs for functional %s", (*f_it).first.c_str());
    // Loop over the functions requested by the form and hook up pointers
    uint ncoeff = (*(*f_it).second).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*(*f_it).second).coefficient_name(i);
      dolfin::info("    Attaching uflsymbol %s", uflsymbol.c_str());
      GenericFunction_ptr function = (*(*system_).bucket()).fetch_uflsymbol(uflsymbol);

      (*(*f_it).second).set_coefficient(uflsymbol, function);
    }
  }

}

// Return a string describing the contents of the function
const std::string FunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

// Return a string describing the contents of functionals_
const std::string FunctionBucket::functionals_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = functionals_.begin(); f_it != functionals_.end(); f_it++ )
  {
    s << indentation << "Functional " << (*f_it).first  << std::endl;
  }
  return s.str();
}

const bool FunctionBucket::include_in_diagnostics() const
{
  dolfin::error("Failed to find virtual function include_in_diagnostics.");
  return false;
}

void FunctionBucket::output()
{
  if (pvdfile_)
  {
    *pvdfile_ << *boost::dynamic_pointer_cast< dolfin::Function >(function());
  }
}

// Empty the function
void FunctionBucket::empty_()
{
  functionals_.clear();
}

