
#include "BoostTypes.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemsWrapper.h"
#include "SpudSystem.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "petscsnes.h"

using namespace buckettools;

// Specific constructor
SpudSolverBucket::SpudSolverBucket(std::string optionpath, System* system) : optionpath_(optionpath), SolverBucket(system)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudSolverBucket::~SpudSolverBucket()
{
  PetscErrorCode perr;

  if(type()=="SNES")
  {
    perr = SNESDestroy(snes_); CHKERRV(perr);
  }

}

// Fill the function using spud and assuming a buckettools schema structure
void SpudSolverBucket::fill()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  base_fill_();

  forms_fill_();

  if (type()=="SNES")
  {
    perr = SNESCreate(PETSC_COMM_WORLD, &snes_); CHKERRV(perr);
    buffer.str(""); buffer << name() << "_";
    perr = SNESSetOptionsPrefix(snes_, buffer.str().c_str()); CHKERRV(perr);

    // we always have at least one ksp
    perr = SNESGetKSP(snes_, &ksp_); CHKERRV(perr);
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver";
    ksp_fill_(buffer.str(), ksp_);
  }
  else if (type()=="Picard")
  {
;
  }
  else
  {
    dolfin::error("Unknown solver type.");
  }

}

void SpudSolverBucket::ksp_fill_(const std::string &optionpath, KSP &ksp)
{
  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  std::string iterative_method;
  buffer.str(""); buffer << optionpath << "/iterative_method/name";
  serr = Spud::get_option(buffer.str(), iterative_method); spud_err(buffer.str(), serr);

  perr = KSPSetType(ksp, iterative_method.c_str()); CHKERRV(perr);

  std::string preconditioner;
  buffer.str(""); buffer << optionpath << "/preconditioner/name";
  serr = Spud::get_option(buffer.str(), preconditioner); spud_err(buffer.str(), serr);

  PC_ptr pc( new PC );
  perr = KSPGetPC(ksp, &(*pc)); CHKERRV(perr);
  perr = PCSetType(*pc, preconditioner.c_str()); CHKERRV(perr);

  pcs_.push_back(pc);

  if (preconditioner=="ksp")
  {

    buffer.str(""); buffer << optionpath << "/preconditioner/linear_solver";
    KSP_ptr subksp;
    perr = PCKSPGetKSP(*pc, &(*subksp)); CHKERRV(perr);
    ksp_fill_(buffer.str(), *subksp);
    subksps_.push_back(subksp);

  }
  else if (preconditioner=="fieldsplit")
  {

    perr = PCFieldSplitSetType(*pc, PC_COMPOSITE_MULTIPLICATIVE); CHKERRV(perr);

    buffer.str(""); buffer << optionpath << "/preconditioner/fieldsplit";
    int nsplits = Spud::option_count(buffer.str());
    for (uint i = 0; i < nsplits; i++)
    {
      buffer.str(""); buffer << optionpath << "/preconditioner/fieldsplit[" << i << "]";
       
      buffer.str(""); buffer << optionpath << "/preconditioner/fieldsplit[" << i << "]/linear_solver";
    }



  }  


}

void SpudSolverBucket::base_fill_()
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;

  // What is this field called?
  buffer.str(""); buffer << optionpath() << "/name";
  serr = Spud::get_option(buffer.str(), name_); spud_err(buffer.str(), serr);

  // Get the field type
  buffer.str(""); buffer << optionpath() << "/type/name";
  serr = Spud::get_option(buffer.str(), type_); spud_err(buffer.str(), serr);

}

void SpudSolverBucket::forms_fill_()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  Spud::OptionError serr;
   
  buffer.str(""); buffer << optionpath() << "/type/form";
  int nforms = Spud::option_count(buffer.str());
  for (uint i = 0; i < nforms; i++)
  {
    buffer.str(""); buffer << optionpath() << "/type/form[" << i << "]";
    std::string formoptionpath = buffer.str();

    std::string formname;
    buffer.str(""); buffer << formoptionpath << "/name";
    serr = Spud::get_option(buffer.str(), formname); spud_err(buffer.str(), serr);

    Form_ptr form = ufc_fetch_form((*system_).name(), name(), type(), formname, (*system_).functionspace());
    register_form(form, formname, formoptionpath);

    // Loop over the functions requested by the form and hook up (potentially currently null) pointers
    uint ncoeff = (*form).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*form).coefficient_name(i);
      GenericFunction_ptr function = (*system_).fetch_uflsymbol(uflsymbol);
      if (!function)
      {
        // the function isn't initialised so either we haven't reached it yet or
        // we got to it but weren't able to initialise it because we didn't have its
        // function space available yet... let's see if it's a function
        std::string functionname = (*system_).fetch_uflname(uflsymbol);
        if (Spud::have_option((*dynamic_cast<SpudSystem*>(system_)).optionpath()+"/coefficient::"+functionname+"/type::Function"))
        {
          // yes, it's a coefficient function... so let's take this opportunity to register
          // its functionspace
          if (!(*system_).contains_coefficientspace(functionname))
          {
            FunctionSpace_ptr coefficientspace;
            coefficientspace = ufc_fetch_coefficientspace((*system_).name(), name(), functionname, (*system_).mesh());
            (*system_).register_coefficientspace(coefficientspace, functionname);
          }
        }
      }

      (*form).set_coefficient(uflsymbol, function);
    }
  }

}

// Register a functional in the function
void SpudSolverBucket::register_form(Form_ptr form, std::string name, std::string optionpath)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it != forms_.end())
  {
    // if it does, issue an error
    dolfin::error("Form named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    forms_[name]            = form;
    form_optionpaths_[name] = optionpath;
  }
}

// Return a string describing the contents of the spudfunction
std::string SpudSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Describe the contents of the form_optionpaths_ map
std::string SpudSolverBucket::forms_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = form_optionpaths_.begin(); s_it != form_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

