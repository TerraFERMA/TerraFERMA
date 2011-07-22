
#include "BoostTypes.h"
#include "SolverBucket.h"
#include "SystemBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
SolverBucket::SolverBucket()
{
  // Do nothing
}

// Specific constructor
SolverBucket::SolverBucket(SystemBucket* system) : system_(system)
{
  // Do nothing
}

// Default destructor
SolverBucket::~SolverBucket()
{
  empty_();
}

// Register a form in the solver
void SolverBucket::register_form(Form_ptr form, std::string name)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it != forms_.end())
  {
    // if it does, issue an error
    dolfin::error("Form named \"%s\" already exists in solver.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    forms_[name] = form;
  }
}

// Return a pointer to a form with the given name
Form_ptr SolverBucket::fetch_form(std::string name)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it == forms_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Form named \"%s\" does not exist in solver.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

Form_it SolverBucket::forms_begin()
{
  return forms_.begin();
}

Form_const_it SolverBucket::forms_begin() const
{
  return forms_.begin();
}

Form_it SolverBucket::forms_end()
{
  return forms_.end();
}

Form_const_it SolverBucket::forms_end() const
{
  return forms_.end();
}

void SolverBucket::attach_form_coeffs()
{
  for (Form_it f_it = forms_begin(); f_it != forms_end(); f_it++)
  {
    // Loop over the functions requested by the form and hook up pointers
    uint ncoeff = (*(*f_it).second).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*(*f_it).second).coefficient_name(i);
      GenericFunction_ptr function = (*system_).fetch_uflsymbol(uflsymbol);

      (*(*f_it).second).set_coefficient(uflsymbol, function);
    }
  }

}

// Return a string describing the contents of the function
const std::string SolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Return a string describing the contents of forms_
const std::string SolverBucket::forms_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = forms_.begin(); f_it != forms_.end(); f_it++ )
  {
    s << indentation << "Form " << (*f_it).first  << std::endl;
  }
  return s.str();
}

// Empty the function
void SolverBucket::empty_()
{
  forms_.clear();
}

