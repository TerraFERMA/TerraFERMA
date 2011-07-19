
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

  PC pc;
  perr = KSPGetPC(ksp, &pc); CHKERRV(perr);
  perr = PCSetType(pc, preconditioner.c_str()); CHKERRV(perr);

  if (preconditioner=="ksp")
  {
    buffer.str(""); buffer << optionpath << "/preconditioner/linear_solver";
    KSP subksp;
    perr = PCKSPGetKSP(pc, &subksp); CHKERRV(perr);
    // recurse!
    ksp_fill_(buffer.str(), subksp);
  }
  else if (preconditioner=="fieldsplit")
  {
    buffer.str(""); buffer << optionpath << "/preconditioner";
    pc_fieldsplit_fill_(buffer.str(), pc);
  }  

}

void SpudSolverBucket::pc_fieldsplit_fill_(const std::string &optionpath, PC &pc)
{

  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  std::string ctype;
  buffer.str(""); buffer << optionpath << "/composite_type/name";
  serr = Spud::get_option(buffer.str(), ctype); spud_err(buffer.str(), serr);
  if (ctype == "additive")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE); CHKERRV(perr);
  }
  else if (ctype == "multiplicative")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE); CHKERRV(perr);
  }
  else if (ctype == "symmetric_multiplicative")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE); CHKERRV(perr);
  }
  else if (ctype == "special")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SPECIAL); CHKERRV(perr);
  }
  else if (ctype == "schur")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR); CHKERRV(perr);
  }
  else
  {
    dolfin::error("Unknown PCCompositeType.");
  }

  if (Spud::have_option(optionpath+"/fieldsplit_by_field"))
  {
    buffer.str(""); buffer << optionpath << "/fieldsplit_by_field/fieldsplit";
    int nsplits = Spud::option_count(buffer.str());
    for (uint i = 0; i < nsplits; i++)
    {
      buffer.str(""); buffer << optionpath << "/fieldsplit_by_field/fieldsplit[" << i << "]";
      pc_fieldsplit_by_field_fill_(buffer.str(), pc);
    }

    KSP *subksps;
    PetscInt nsubksps;
    perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); CHKERRV(perr); 

    assert(nsubksps==nsplits);

    for (uint i = 0; i < nsplits; i++)
    {
      buffer.str(""); buffer << optionpath << "/fieldsplit_by_field/fieldsplit[" << i << "]/linear_solver";
      // recurse!
      ksp_fill_(buffer.str(), subksps[i]);
    }

  }
  else if (Spud::have_option(optionpath+"/fieldsplit_by_region"))
  {
    buffer.str(""); buffer << optionpath << "/fieldsplit_by_region/fieldsplit";
    int nsplits = Spud::option_count(buffer.str());
    for (uint i = 0; i < nsplits; i++)
    {
      buffer.str(""); buffer << optionpath << "/fieldsplit_by_region/fieldsplit[" << i << "]";
      pc_fieldsplit_by_region_fill_(buffer.str(), pc);
    }
  }
  else
  {
    dolfin::error("Unknown way of specifying fieldsplit");
  }

  

}

void SpudSolverBucket::pc_fieldsplit_by_field_fill_(const std::string &optionpath, PC &pc)
{

  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  std::vector< uint > indices_vector;

  buffer.str(""); buffer << optionpath << "/field";
  int nfields = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfields; i++)
  {
    buffer.str(""); buffer << optionpath << "/field[" << i << "]";
    std::string fieldname;
    buffer.str(""); buffer << optionpath << "/field[" << i << "]/name";
    serr = Spud::get_option(buffer.str(), fieldname); spud_err(buffer.str(), serr);
    
    int fieldindex = (*(*system_).fetch_field(fieldname)).index();

    buffer.str(""); buffer << optionpath << "/field[" << i << "]/components";
    if (Spud::have_option(buffer.str()))
    {
      std::vector< int > components;
      serr = Spud::get_option(buffer.str(), components); spud_err(buffer.str(), serr);
      
      for (std::vector<int>::const_iterator comp = components.begin(); comp != components.end(); comp++)
      {
        assert(*comp < (*(*(*system_).functionspace())[fieldindex]).element().num_sub_elements());

        boost::unordered_set<uint> dof_set = (*(*(*(*system_).functionspace())[fieldindex])[*comp]).dofmap().dofs();
        std::pair<uint, uint> ownership_range = (*(*system_).functionspace()).dofmap().ownership_range();

        for (boost::unordered_set<uint>::const_iterator dof_it = dof_set.begin(); dof_it != dof_set.end(); dof_it++)
        {
          if ((*dof_it >= ownership_range.first) && (*dof_it < ownership_range.second))
          {
            indices_vector.push_back(*dof_it);
          }
        }
      }
    }
    else
    {
      boost::unordered_set<uint> dof_set = (*(*(*system_).functionspace())[fieldindex]).dofmap().dofs();
      std::pair<uint, uint> ownership_range = (*(*system_).functionspace()).dofmap().ownership_range();

      for (boost::unordered_set<uint>::const_iterator dof_it = dof_set.begin(); dof_it != dof_set.end(); dof_it++)
      {
        if ((*dof_it >= ownership_range.first) && (*dof_it < ownership_range.second))
        {
          indices_vector.push_back(*dof_it);
        }
      }
    }
  }

  std::sort(indices_vector.begin(), indices_vector.end());

  PetscInt n=indices_vector.size();
  PetscInt *indices;

  PetscMalloc(n*sizeof(PetscInt), &indices);
  uint ind = 0;
  for (std::vector<uint>::const_iterator ind_it = indices_vector.begin(); ind_it != indices_vector.end(); ind_it++)
  {
    indices[ind] = *ind_it;
    ind++;
  }
  assert(ind==n);
  IS is;
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, &is); CHKERRV(perr);
  perr = ISView(is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);
  perr = PCFieldSplitSetIS(pc, is); CHKERRV(perr);

  PetscFree(indices);
   
}

void SpudSolverBucket::pc_fieldsplit_by_region_fill_(const std::string &optionpath, PC &pc)
{

  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

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

