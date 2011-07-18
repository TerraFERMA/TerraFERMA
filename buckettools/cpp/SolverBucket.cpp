
#include "BoostTypes.h"
#include "SolverBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
SolverBucket::SolverBucket()
{
  // Do nothing
}

// Specific constructor
SolverBucket::SolverBucket(System* system) : system_(system)
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

// Return a string describing the contents of the function
std::string SolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Return a string describing the contents of forms_
std::string SolverBucket::forms_str(int indent) const
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

