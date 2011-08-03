
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

const FunctionBucket_ptr SystemBucket::fetch_field(std::string name) const
{
  // First check if a field with this name already exists
  FunctionBucket_const_it f_it = fields_.find(name);
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

FunctionBucket_it SystemBucket::fields_begin()
{
  return fields_.begin();
}

FunctionBucket_const_it SystemBucket::fields_begin() const
{
  return fields_.begin();
}

FunctionBucket_it SystemBucket::fields_end()
{
  return fields_.end();
}

FunctionBucket_const_it SystemBucket::fields_end() const
{
  return fields_.end();
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

FunctionBucket_it SystemBucket::coeffs_begin()
{
  return coeffs_.begin();
}

FunctionBucket_const_it SystemBucket::coeffs_begin() const
{
  return coeffs_.begin();
}

FunctionBucket_it SystemBucket::coeffs_end()
{
  return coeffs_.end();
}

FunctionBucket_const_it SystemBucket::coeffs_end() const
{
  return coeffs_.end();
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
    // for the moment assume that the order they are added is the order you want them solved
    // (this is easily extendible to a user defined order)
    orderedsolvers_[(int) solvers_.size()] = solver;
  }
}

SolverBucket_it SystemBucket::solvers_begin()
{
  return solvers_.begin();
}

SolverBucket_const_it SystemBucket::solvers_begin() const
{
  return solvers_.begin();
}

SolverBucket_it SystemBucket::solvers_end()
{
  return solvers_.end();
}

SolverBucket_const_it SystemBucket::solvers_end() const
{
  return solvers_.end();
}

int_SolverBucket_it SystemBucket::orderedsolvers_begin()
{
  return orderedsolvers_.begin();
}

int_SolverBucket_const_it SystemBucket::orderedsolvers_begin() const
{
  return orderedsolvers_.begin();
}

int_SolverBucket_it SystemBucket::orderedsolvers_end()
{
  return orderedsolvers_.end();
}

int_SolverBucket_const_it SystemBucket::orderedsolvers_end() const
{
  return orderedsolvers_.end();
}

void SystemBucket::solve()
{
  for (int_SolverBucket_const_it s_it = orderedsolvers_begin(); s_it != orderedsolvers_end(); s_it++)
  {
    (*(*s_it).second).solve();
  }
}

void SystemBucket::output()
{
  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end(); f_it++)
  {
    (*(*f_it).second).output();
  }
}

void SystemBucket::collect_bcs_()
{
  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end(); f_it++)
  {
    for (BoundaryCondition_const_it b_it = (*(*f_it).second).bcs_begin(); b_it != (*(*f_it).second).bcs_end(); b_it++)
    {
      bcs_.push_back((*b_it).second);
    }
  }
}

std::vector<BoundaryCondition_ptr>::iterator SystemBucket::bcs_begin()
{
  return bcs_.begin();
}

std::vector<BoundaryCondition_ptr>::const_iterator SystemBucket::bcs_begin() const
{
  return bcs_.begin();
}

std::vector<BoundaryCondition_ptr>::iterator SystemBucket::bcs_end()
{
  return bcs_.end();
}

std::vector<BoundaryCondition_ptr>::const_iterator SystemBucket::bcs_end() const
{
  return bcs_.end();
}

void SystemBucket::attach_and_initialize()
{
  // attach the coefficients to the functionals and forms
  // now that we have them all registered
  dolfin::info("Attaching coeffs for system %s", name().c_str());
  attach_all_coeffs_();

  for (SolverBucket_it s_it = solvers_begin(); s_it != solvers_end(); s_it++)
  {
    (*(*s_it).second).initialize_matrices();
  }

}

void SystemBucket::attach_all_coeffs_()
{
  attach_function_coeffs_(fields_begin(), fields_end());
  attach_function_coeffs_(coeffs_begin(), coeffs_end());
  attach_solver_coeffs_(solvers_begin(), solvers_end());
}

void SystemBucket::attach_function_coeffs_(FunctionBucket_it f_begin, FunctionBucket_it f_end)
{
  for (FunctionBucket_it f_it = f_begin; f_it != f_end; f_it++)
  {
    (*(*f_it).second).attach_functional_coeffs();
  }
}

void SystemBucket::attach_solver_coeffs_(SolverBucket_it s_begin, SolverBucket_it s_end)
{

  for (SolverBucket_it s_it = s_begin; s_it != s_end; s_it++)
  {
    (*(*s_it).second).attach_form_coeffs();
  }
}

// Return a string describing the contents of the system
const std::string SystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << std::endl;
  indent++;
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  return s.str();
}

// Return a string describing the contents of fields_
const std::string SystemBucket::fields_str(int indent) const
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
const std::string SystemBucket::coeffs_str(int indent) const
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
const std::string SystemBucket::solvers_str(int indent) const
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

