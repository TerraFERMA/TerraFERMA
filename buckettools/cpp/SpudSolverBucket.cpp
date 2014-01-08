// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#include "BoostTypes.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemSolversWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "petscsnes.h"
#include "ConvergenceFile.h"
#include "KSPConvergenceFile.h"
#include "PythonExpression.h"
#include <boost/algorithm/string/predicate.hpp>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudSolverBucket::SpudSolverBucket(const std::string &optionpath, 
                                            SystemBucket* system) : 
                                            optionpath_(optionpath), 
                                            SolverBucket(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudSolverBucket::~SpudSolverBucket()
{
  empty_();                                                          // empty the data in the derived class
}

//*******************************************************************|************************************************************//
// fill the solver bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill()
{

  fill_base_();                                                      // fill the base solver data: type, name etc.

  fill_forms_();                                                     // fill the forms data

}

//*******************************************************************|************************************************************//
// initialize the actual solvers in the solver bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSolverBucket::initialize()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  initialize_tensors_();                                             // set up the tensor structures

  std::stringstream prefix;                                          // prefix buffer
  prefix.str(""); prefix << (*system_).name() << "_" << name() << "_";

  if (type()=="SNES")                                                // if this is a snes solver.  FIXME: switch to enum check
  {

    perr = SNESCreate(PETSC_COMM_WORLD, &snes_); CHKERRV(perr);      // create the petsc snes object

    perr = SNESSetOptionsPrefix(snes_, prefix.str().c_str());        // set its petsc options name prefix to SystemName_SolverName
    CHKERRV(perr);

    ctx_.solver = this;                                              // the snes context just needs this class... neat, huh?

    perr = SNESSetFunction(snes_, *(*res_).vec(),                    // set the snes function to use the newly allocated residual vector
                                    FormFunction, (void *) &ctx_); 
    CHKERRV(perr);

    if (bilinearpc_)                                                 // if we have a pc form
    {
      assert(matrixpc_);
      perr = SNESSetJacobian(snes_, *(*matrix_).mat(),               // set the snes jacobian to have two matrices
                  *(*matrixpc_).mat(), FormJacobian, (void *) &ctx_); 
      CHKERRV(perr);
    }
    else                                                             // otherwise
    {
      perr = SNESSetJacobian(snes_, *(*matrix_).mat(),               // set the snes jacobian to have the same matrix twice
                    *(*matrix_).mat(), FormJacobian, (void *) &ctx_); 
      CHKERRV(perr);
    }

    std::string snestype;
    buffer.str(""); buffer << optionpath() << "/type/snes_type/name";
    serr = Spud::get_option(buffer.str(), snestype);                 // set the snes type... ls is most common
    spud_err(buffer.str(), serr);

    if(snestype=="vi")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #if PETSC_VERSION_MINOR > 3
      perr = SNESSetType(snes_, SNESVINEWTONRSLS); CHKERRV(perr); 
      #elif PETSC_VERSION_MINOR == 3
      perr = SNESSetType(snes_, SNESVIRS); CHKERRV(perr); 
      #else
      perr = SNESSetType(snes_, snestype.c_str()); CHKERRV(perr); 
      #endif
      perr = SNESSetFromOptions(snes_); CHKERRV(perr);               // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)

      fill_constraints_();
      #else
      dolfin::error("Cannot set snes vi with PETSc < 3.2.");
      #endif
    }
    else if(snestype=="ls")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
      perr = SNESSetType(snes_, SNESNEWTONLS); CHKERRV(perr);
      #else
      perr = SNESSetType(snes_, snestype.c_str()); CHKERRV(perr);
      #endif 
      perr = SNESSetFromOptions(snes_); CHKERRV(perr);               // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)

      std::string lstype;
      buffer.str(""); buffer << optionpath() << "/type/snes_type::ls/ls_type/name";
      serr = Spud::get_option(buffer.str(), lstype);                // set the snes type... cubic is the most common
      spud_err(buffer.str(), serr);                                 // FIXME: once dev is released (PETSc>3.2) the
                                                                    // schema should be refactored to offer the petsc options
                                                                    // directly rather than the mapping onto options done here
      if (lstype=="cubic")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); CHKERRV(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); CHKERRV(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "bt"); CHKERRV(perr);
        perr = SNESLineSearchSetOrder(linesearch, 3); CHKERRV(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchCubic, PETSC_NULL); CHKERRV(perr); 
        #endif
      }
      else if (lstype=="quadratic")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); CHKERRV(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); CHKERRV(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "bt"); CHKERRV(perr);
        perr = SNESLineSearchSetOrder(linesearch, 2); CHKERRV(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchQuadratic, PETSC_NULL); CHKERRV(perr); 
        #endif
      }
      else if (lstype=="basic")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); CHKERRV(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); CHKERRV(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "basic"); CHKERRV(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchNo, PETSC_NULL); CHKERRV(perr); 
        #endif
      }
      else if (lstype=="basicnonorms")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        dolfin::error("No equivalent snes ls type to basicnonorms in petsc>3.2.");
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchNoNorms, PETSC_NULL); CHKERRV(perr); 
        #endif
      }
      else
      {
        dolfin::error("Unknown snes ls type.");
      }

                                                                     // FIXME: realign with available PETSc options once petsc-dev
                                                                     // is released (PETSc>3.2)
      buffer.str(""); buffer << optionpath() << "/type/snes_type::ls/alpha";
      double alpha;
      serr = Spud::get_option(buffer.str(), alpha, 1.e-4);
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << optionpath() << "/type/snes_type::ls/max_step";
      double maxstep;
      serr = Spud::get_option(buffer.str(), maxstep, 1.e8);
      spud_err(buffer.str(), serr);
       
      buffer.str(""); buffer << optionpath() << "/type/snes_type::ls/min_lambda";
      if (Spud::have_option(buffer.str()))
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        dolfin::error("Cannot set snes ls min_lambda with PETSc > 3.2 - options not aligned.");
        #endif
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 2
        dolfin::error("Cannot set snes ls min_lambda with PETSc < 3.2.");
        #endif
      }
      double minlambda;
      serr = Spud::get_option(buffer.str(), minlambda, 1.e-12);
      spud_err(buffer.str(), serr);

      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #if PETSC_VERSION_MINOR > 2
      SNESLineSearch linesearch;
      #if PETSC_VERSION_MINOR > 3
      perr = SNESGetLineSearch(snes_, &linesearch); CHKERRV(perr);
      #else
      perr = SNESGetSNESLineSearch(snes_, &linesearch); CHKERRV(perr);
      #endif
      perr = SNESLineSearchBTSetAlpha(linesearch, alpha); CHKERRV(perr);// FIXME: assumes using bt
      perr = SNESLineSearchSetTolerances(linesearch, PETSC_DEFAULT, maxstep, PETSC_DEFAULT,
                                        PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
      CHKERRV(perr);
      #else
      perr = SNESLineSearchSetParams(snes_, alpha, maxstep, minlambda); CHKERRV(perr);
      #endif
      #else
      perr = SNESLineSearchSetParams(snes_, alpha, maxstep); CHKERRV(perr);
      #endif
       
    }
    else
    {
      perr = SNESSetType(snes_, snestype.c_str()); CHKERRV(perr); 
      perr = SNESSetFromOptions(snes_); CHKERRV(perr);               // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)
    }

    buffer.str(""); buffer << optionpath() 
                                        << "/type/monitors/residual";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorDefault,               // set a snes residual monitor
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath() 
                                  << "/type/monitors/solution_graph";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorSolution,              // set a snes solution monitor (graph)
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
    {
      snesmctx_.solver = this;
      if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
      {
        buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                               << (*system()).name() << "_" 
                               << name() << "_snes.conv";
        convfile_.reset( new ConvergenceFile(buffer.str(),           // allocate the file but don't write the header yet as the
                                      (*system()).name(), name()) ); // bucket isn't complete
      }
      perr = SNESMonitorSet(snes_, SNESCustomMonitor,                // set a custom snes monitor
                                            &snesmctx_, PETSC_NULL); 
      CHKERRV(perr);
    }

    perr = SNESSetTolerances(snes_, atol_, rtol_, stol_, maxits_,    // from the data we collected in the base data fill set-up the
                                                            maxfes_);// snes tolerances
    CHKERRV(perr);

    perr = SNESGetKSP(snes_, &ksp_); CHKERRV(perr);                  // we always have at least one ksp so use the solverbucket ksp
                                                                     // to start setting up the ksp inside the snes
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // the ksp solver path
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // can then be used to fill the ksp data

    buffer.str(""); buffer << optionpath() << "/type/monitors/view_snes";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESView(snes_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);// turn on snesview so we get some debugging info
    }

  }
  else if (type()=="Picard")                                         // if this is a picard solver
  {

    perr = KSPCreate(PETSC_COMM_WORLD, &ksp_); CHKERRV(perr);        // create a ksp object from the variable in the solverbucket

    if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
    {
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_" 
                             << name() << "_picard.conv";
      convfile_.reset( new ConvergenceFile(buffer.str(),             // allocate the file but don't write the header yet as the
                                    (*system()).name(), name()) );   // bucket isn't complete
    }

    if (bilinearpc_)
    {                                                                // if there's a pc associated
      perr = KSPSetOperators(ksp_, *(*matrix_).mat(),                // set the ksp operators with two matrices
                                   *(*matrixpc_).mat(), 
                                   SAME_NONZERO_PATTERN); 
      CHKERRV(perr);
    }
    else
    {
      perr = KSPSetOperators(ksp_, *(*matrix_).mat(),                // set the ksp operators with the same matrices
                                   *(*matrix_).mat(), 
                                   SAME_NONZERO_PATTERN); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // figure out the linear solver optionspath
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // fill the ksp data

    buffer.str(""); buffer << optionpath() << "/type/linear_solver/monitors/view_ksp";
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPView(ksp_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr); // turn on kspview so we get some debugging info
    }

  }
  else                                                               // don't know how we got here
  {
    dolfin::error("Unknown solver type.");
  }

}

//*******************************************************************|************************************************************//
// make a partial copy of the provided solver bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SpudSolverBucket::copy_diagnostics(SolverBucket_ptr &solver, SystemBucket_ptr &system) const
{

  if(!solver)
  {
    solver.reset( new SpudSolverBucket(optionpath_, &(*system)) );
  }

  SolverBucket::copy_diagnostics(solver, system);

  (*boost::dynamic_pointer_cast< SpudSolverBucket >(solver)).form_optionpaths_ = form_optionpaths_;

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps (with an optionpath as well)
//*******************************************************************|************************************************************//
void SpudSolverBucket::register_form(Form_ptr form, 
                                      const std::string &name, 
                                      const std::string &optionpath)
{
  Form_it f_it = forms_.find(name);                                  // check if name exists
  if (f_it != forms_.end())
  {
    dolfin::error("Form named \"%s\" already exists in function.",   // if it does, issue an error
                                                    name.c_str());
  }
  else
  {
    forms_[name]            = form;                                  // if not, insert form pointer into data map
    form_optionpaths_[name] = optionpath;                            // and do the same for its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the solver bucket
//*******************************************************************|************************************************************//
const std::string SpudSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << " (" << 
                                    optionpath() << ")" << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the forms in the solver bucket
//*******************************************************************|************************************************************//
const std::string SpudSolverBucket::forms_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = form_optionpaths_.begin(); 
                            s_it != form_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << 
                                (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// fill the solver bucket base data assuming the buckettools schema (common for all solver types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_base_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // solver name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/name";            // solver type (as string)
  serr = Spud::get_option(buffer.str(), type_);                      // FIXME: add conversion to enum here
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/relative_error";  // relative nonlinear error
  serr = Spud::get_option(buffer.str(), rtol_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/absolute_error";  // absolute nonlinear error
  serr = Spud::get_option(buffer.str(), atol_, 1.e-50); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/solution_error";  // nonlinear solution error (only applies to snes solver types)
  serr = Spud::get_option(buffer.str(), stol_, 1.e-8); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/max_iterations";  // maximum number of nonlinear iterations
  serr = Spud::get_option(buffer.str(), maxits_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/min_iterations";  // minimum number of nonlinear iterations (only applies to
  serr = Spud::get_option(buffer.str(), minits_, 0);                 // picard solver types)
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() <<                          // maximum number of residual evaluations (only applies to snes
                                  "/type/max_function_evaluations";  // solver types)
  serr = Spud::get_option(buffer.str(), maxfes_, 10000); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << 
                                "/type/ignore_all_solver_failures";
  ignore_failures_ = Spud::have_option(buffer.str());

  iteration_count_.reset( new int );
  *iteration_count_ = 0;

  copy_ = false;

  buffer.str(""); buffer << optionpath() << 
                                "/type/monitors/norms";
  monitornorms_ = Spud::have_option(buffer.str());

}

//*******************************************************************|************************************************************//
// fill the forms in solver bucket assuming the buckettools schema (common for all solver types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_forms_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
   
  buffer.str(""); buffer << optionpath() << "/type";
  fill_subforms_(buffer.str());
  std::stringstream prefix;                                          // prefix buffer
  prefix.str(""); prefix << (*system_).name() << "_" << name() << "_";
  fill_solverforms_(buffer.str(), prefix.str());

                                                                     // depending on the type of solver we assign certain forms to
                                                                     // hardcoded names in the solver bucket (basically it depends
                                                                     // which are linear or bilinear)
  ident_zeros_ = false;
  ident_zeros_pc_ = false;

  if (type()=="SNES")                                                // snes solver type...
  {
    linear_      = fetch_form("Residual");
    bilinear_    = fetch_form("Jacobian");
    buffer.str(""); buffer << optionpath() << "/type::SNES/form::Jacobian/ident_zeros";
    ident_zeros_ = Spud::have_option(buffer.str());


    if (contains_form("JacobianPC"))                                 // is there a pc form?
    {
      bilinearpc_ = fetch_form("JacobianPC");                        // yes but
      buffer.str(""); buffer << optionpath() << "/type::SNES/form::JacobianPC/ident_zeros";
      ident_zeros_pc_ = Spud::have_option(buffer.str());
    }                                                                // otherwise bilinearpc_ is null (indicates self pcing)
    assert(!residual_);                                              // residual_ is always a null pointer for snes

  }
  else if (type()=="Picard")                                         // picard solver type...
  {
    linear_      = fetch_form("Linear");
    bilinear_    = fetch_form("Bilinear");
    buffer.str(""); buffer << optionpath() << "/type::Picard/form::Bilinear/ident_zeros";
    ident_zeros_ = Spud::have_option(buffer.str());


    if (contains_form("BilinearPC"))                                 // is there a pc form?
    {
      bilinearpc_ = fetch_form("BilinearPC");                        // yes but
      buffer.str(""); buffer << optionpath() << "/type::Picard/form::BilinearPC/ident_zeros";
      ident_zeros_pc_ = Spud::have_option(buffer.str());
    }                                                                // otherwise bilinearpc_ is null (indicates self pcing)
    residual_   = fetch_form("Residual");

  }
  else                                                               // unknown solver type
  {
    dolfin::error("Unknown solver type.");
  }

  for (Form_const_it f_it = forms_begin(); f_it != forms_end(); f_it++)
  {
    if(boost::algorithm::ends_with((*f_it).first, "SchurPC"))
    {
      register_solverform((*f_it).second, (*f_it).first);
      solverident_zeros_[(*f_it).first] = Spud::have_option(form_optionpaths_[(*f_it).first]+"/ident_zeros");
    }
  }

}

//*******************************************************************|************************************************************//
// fill a form in solver bucket assuming the buckettools schema (common for all solver types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_subforms_(const std::string &optionpath, 
                                      const std::string &prefix)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
   
  buffer.str(""); buffer << optionpath << "/form";           
  int nforms = Spud::option_count(buffer.str());                     // find out how many forms we have
  for (uint i = 0; i < nforms; i++)                                  // loop over them 
  {
    buffer.str(""); buffer << optionpath << "/form[" << i << "]";
    std::string formoptionpath = buffer.str();

    std::string formname;
    buffer.str(""); buffer << formoptionpath << "/name";             // get the form name
    serr = Spud::get_option(buffer.str(), formname); 
    spud_err(buffer.str(), serr);

    Form_ptr form = ufc_fetch_form((*system_).name(), name(),        // get the form from the ufc using the system functionspace
                                        type(), 
                                        prefix+formname, 
                                        (*system_).functionspace());
    (*form).set_cell_domains((*system_).celldomains());
    (*form).set_interior_facet_domains((*system_).facetdomains());
    (*form).set_exterior_facet_domains((*system_).facetdomains());
    register_form(form, prefix+formname, formoptionpath);

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // form because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
    uint ncoeff = (*form).num_coefficients();                        // how many coefficients does this form require?
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*form).coefficient_name(i);           // what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
      if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))  // a base ufl symbol was only inserted into the parent bucket's
      {                                                              // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
        std::string baseuflsymbol =                                  // what is the base ufl symbol?
              (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol); // have we already registered a functionspace for this base ufl
                                                                     // symbol?
        if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
        {                                                            // no...
          FunctionSpace_ptr coefficientspace;
          coefficientspace = ufc_fetch_coefficientspace_from_solver( // take a pointer to the functionspace from the ufc
                                        (*system_).name(), name(), 
                                        baseuflsymbol, 
                                        (*system_).mesh());
          (*(*system_).bucket()).register_coefficientspace(          // and register it in the parent bucket's map
                                        coefficientspace, 
                                        baseuflsymbol);
        }
      }

    }

  }

}

//*******************************************************************|************************************************************//
// recursively fill any forms beneath the linear_solver
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_solverforms_(const std::string &optionpath,
                                         const std::string &prefix)
{
  std::stringstream buffer, fsbuffer, lsbuffer, pcbuffer;            // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
   
  buffer.str(""); buffer << optionpath <<
                 "/composite_type::schur/schur_preconditioner::user";
  if (Spud::have_option(buffer.str()))
  {
    fill_subforms_(buffer.str(), prefix);
  }

  fsbuffer.str(""); fsbuffer << optionpath << "/fieldsplit";
  lsbuffer.str(""); lsbuffer << optionpath <<
                                     "/linear_solver/preconditioner";
  pcbuffer.str(""); pcbuffer << optionpath << "preconditioner";
  if (Spud::have_option(fsbuffer.str()))
  {
    for (uint i = 0; i < Spud::option_count(fsbuffer.str()); i++)
    {
      std::string fsname;
      buffer.str(""); buffer << fsbuffer.str() << "[" << i << "]/name";
      serr = Spud::get_option(buffer.str(), fsname);
      spud_err(buffer.str(), serr);
      buffer.str(""); buffer << fsbuffer.str() << "[" << i << "]/linear_solver/preconditioner";
      fill_solverforms_(buffer.str(), prefix+fsname+"_");
    }
  }
  else if (Spud::have_option(lsbuffer.str()))
  {
    fill_solverforms_(lsbuffer.str(), prefix);
  }
  else if (Spud::have_option(pcbuffer.str()))
  {
    fill_solverforms_(pcbuffer.str(), prefix);
  }

}

//*******************************************************************|************************************************************//
// initialize the tensor structures in solver bucket assuming the buckettools schema (common for all solver types)
// (must be called after fill_forms_)
//*******************************************************************|************************************************************//
void SpudSolverBucket::initialize_tensors_()
{
  std::stringstream buffer;                                          // optionpath buffer
  dolfin::AssemblerBase assembler;
  assembler.finalize_tensor = false;
  assembler.keep_diagonal = true;
   
  work_.reset( new dolfin::PETScVector(*boost::dynamic_pointer_cast<dolfin::PETScVector>((*(*system_).function()).vector())) ); 
  (*work_).zero();

  matrix_.reset(new dolfin::PETScMatrix);                            // allocate the matrix
  assembler.init_global_tensor(*matrix_, *bilinear_);

  if(bilinearpc_)                                                    // do we have a pc form?
  {
    matrixpc_.reset(new dolfin::PETScMatrix);                        // allocate the matrix
    assembler.init_global_tensor(*matrixpc_, *bilinearpc_);
  }
     
  rhs_.reset(new dolfin::PETScVector);                               // allocate the rhs
  assembler.init_global_tensor(*rhs_, *linear_);

  rhsbc_.reset(new dolfin::PETScVector);                             // allocate the rhs
  assembler.init_global_tensor(*rhsbc_, *linear_);

  if(residual_)                                                      // do we have a residual_ form?
  {                                                                  // yes...
    res_.reset(new dolfin::PETScVector);                             // allocate the residual
    assembler.init_global_tensor(*res_, *residual_);
  }
  else
  {
    res_.reset( new dolfin::PETScVector(*work_) );                   // but we still want to initialize the residual vector
    (*res_).zero();
  }

  for (Form_const_it f_it = solverforms_begin(); 
                     f_it != solverforms_end(); f_it++)
  {
    PETScMatrix_ptr solvermatrix;
    solvermatrix.reset(new dolfin::PETScMatrix);
    assembler.init_global_tensor(*solvermatrix, 
                                 *(*f_it).second);
    solvermatrices_[(*f_it).first] = solvermatrix;
  }

  assemble_linearforms();                                            // preassemble
  assemble_bilinearforms();

}

//*******************************************************************|************************************************************//
// fill a ksp object from the options tree (not necessarily the main solver bucket ksp_ object as this routine may be called
// recursively for ksp and fieldsplit pc types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                 const std::string prefix,
                                 const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  perr = KSPSetOptionsPrefix(ksp, prefix.c_str());                   // set the ksp options prefix
  CHKERRV(perr);

  std::string iterative_method;
  buffer.str(""); buffer << optionpath << "/iterative_method/name";  // iterative method (gmres, fgmres, cg etc.)
  serr = Spud::get_option(buffer.str(), iterative_method); 
  spud_err(buffer.str(), serr);

  perr = KSPSetType(ksp, iterative_method.c_str()); CHKERRV(perr);   // set the ksp type to the iterative method

  perr = KSPSetFromOptions(ksp);                                     // do this now so that options will be overwritten by options
                                                                     // file

  PC pc;
  perr = KSPGetPC(ksp, &pc); CHKERRV(perr);                          // get the pc from the ksp

  fill_pc_(optionpath, pc,
           prefix,
           parent_indices);

  if(iterative_method != "preonly")                                  // tolerances (and monitors) only apply to iterative methods
  {
    PetscReal rtol;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/relative_error";
    serr = Spud::get_option(buffer.str(), rtol);                     // relative error
    spud_err(buffer.str(), serr);

    PetscReal atol;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/absolute_error";
    serr = Spud::get_option(buffer.str(), atol, 1.e-50);             // absolute error
    spud_err(buffer.str(), serr);

    PetscReal dtol;
    buffer.str(""); buffer << optionpath << 
                                "/iterative_method/divergence_error";
    serr = Spud::get_option(buffer.str(), dtol, 10000.0);            // divergence error (tolerance?) - preconditioned solution must
    spud_err(buffer.str(), serr);                                    // diverge this much to be considered divergent

    PetscInt maxits;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/max_iterations";// maximum number of linear solver iterations
    serr = Spud::get_option(buffer.str(), maxits); 
    spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << 
                              "/iterative_method/zero_initial_guess";// starting the solve from zero?
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE); 
      CHKERRV(perr);
    }
    else
    {
      perr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << 
                "/iterative_method/monitors/preconditioned_residual";// monitor the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorDefault, 
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << 
                          "/iterative_method/monitors/true_residual";// monitor the true residual (more expensive)
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, 
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << 
          "/iterative_method/monitors/preconditioned_residual_graph";// plot a graph of the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
      perr = KSPMonitorSet(ksp, KSPMonitorLGResidualNorm, 
                                             PETSC_NULL, PETSC_NULL); 
      #else
      perr = KSPMonitorSet(ksp, KSPMonitorLG, 
                                             PETSC_NULL, PETSC_NULL); 
      #endif
      CHKERRV(perr);
    }

    if (Spud::have_option(optionpath+"/iterative_method/monitors/convergence_file"))
    {
      kspmctx_.solver = this;
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_" 
                             << name() << "_ksp.conv";
      kspconvfile_.reset( new KSPConvergenceFile(buffer.str(),       // allocate the file but don't write the header yet as the
                                    (*system()).name(), name()) );   // bucket isn't complete
      perr = KSPMonitorSet(ksp, KSPCustomMonitor, 
                                             &kspmctx_, PETSC_NULL); 
      CHKERRV(perr);
    }

    if (Spud::have_option(optionpath+"/iterative_method/monitors/test_null_space"))
    {
      kspmctx_.solver = this;
      perr = KSPMonitorSet(ksp, KSPNullSpaceMonitor, 
                                             &kspmctx_, PETSC_NULL); 
      CHKERRV(perr);
    }

    perr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  }

  buffer.str(""); buffer << optionpath << "/remove_null_space";      // removing a (or multiple) null space(s)
  if (Spud::have_option(buffer.str()))
  {
    MatNullSpace SP;                                                 // create a set of nullspaces in a null space object
    fill_nullspace_(buffer.str(), SP, parent_indices);

    perr = KSPSetNullSpace(ksp, SP); CHKERRV(perr);                  // attach it to the ksp
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatNullSpaceDestroy(&SP); CHKERRV(perr);                  // destroy the null space object, necessary?
    #else
    perr = MatNullSpaceDestroy(SP); CHKERRV(perr);                   // destroy the null space object, necessary?
    #endif
  }

}

//*******************************************************************|************************************************************//
// fill a pc object from the options tree (this routine may be called recursively for ksp, bjacobi, asm and fieldsplit pc types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_pc_(const std::string &optionpath, PC &pc, 
                                const std::string prefix,
                                const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  perr = PCSetOptionsPrefix(pc, prefix.c_str());                     // set the pc options prefix
  CHKERRV(perr);

  std::string preconditioner;
  buffer.str(""); buffer << optionpath << "/preconditioner/name";    // preconditioner type
  serr = Spud::get_option(buffer.str(), preconditioner);
  spud_err(buffer.str(), serr);

  perr = PCSetType(pc, preconditioner.c_str()); CHKERRV(perr);       // set its type (read from options earlier)

  perr = PCSetFromOptions(pc); CHKERRV(perr);                        // do this now so they can be overwritten

  if (preconditioner=="ksp")                                         // if the pc is itself a ksp
  {
    buffer.str(""); buffer << optionpath << 
                                    "/preconditioner/linear_solver";
    KSP subksp;                                                      // create a subksp from this pc
    perr = PCKSPGetKSP(pc, &subksp); CHKERRV(perr);
    fill_ksp_(buffer.str(), subksp, prefix, parent_indices);// recursively fill the ksp data (i.e. go back to this routine)
  }
  else if (preconditioner=="lsc")                                    // would be nice to put subksp options for lsc here
  {
    perr = PCSetUp(pc); CHKERRV(perr);                               // this is just necessary in case we view it later
  }
  else if (preconditioner=="fieldsplit")                             // if the pc is a fieldsplit
  {
    buffer.str(""); buffer << optionpath << "/preconditioner";
    fill_pc_fieldsplit_(buffer.str(), pc, prefix, parent_indices);   // fill the fieldsplit data (will end up back here again
                                                                     // eventually)
  }
  else if (preconditioner=="lu")                                     // if the pc is direct
  {
    std::string factorization_package;                               // we get to choose a factorization package
    buffer.str(""); buffer << optionpath << 
                        "/preconditioner/factorization_package/name";
    serr = Spud::get_option(buffer.str(), factorization_package); 
    spud_err(buffer.str(), serr);

    perr = PCFactorSetMatSolverPackage(pc, 
                                      factorization_package.c_str()); 
    CHKERRV(perr);

  }
  else if (preconditioner=="hypre")
  {
    std::string hypre_type;                                          // we get to choose the hypre type
    buffer.str(""); buffer << optionpath << 
                        "/preconditioner/hypre_type/name";
    serr = Spud::get_option(buffer.str(), hypre_type); 
    spud_err(buffer.str(), serr);

#ifdef PETSC_HAVE_HYPRE
    perr = PCHYPRESetType(pc, hypre_type.c_str()); CHKERRV(perr);
#else
    dolfin::error("Must compile petsc with hypre to use it.");
#endif
  }
  else if ((preconditioner=="bjacobi")||(preconditioner=="asm"))
  {
    perr = PCSetUp(pc); CHKERRV(perr);                               // call this before subpc can be retrieved

    KSP *subksp;
    if(preconditioner=="bjacobi")
    {
      perr = PCBJacobiGetSubKSP(pc, PETSC_NULL, PETSC_NULL, &subksp);
      CHKERRV(perr);
    }
    else if (preconditioner=="asm")
    {
      perr = PCASMGetSubKSP(pc, PETSC_NULL, PETSC_NULL, &subksp);
      CHKERRV(perr);
    }
    else
    {
      dolfin::error("Unknown preconditioner.");
    }

    PC subpc;
    perr = KSPGetPC(*subksp, &subpc); CHKERRV(perr);                  // get the sub pc from the sub ksp
    fill_pc_(optionpath+"/preconditioner", subpc,
             prefix, 
             parent_indices);

  }

  buffer.str(""); buffer << optionpath 
                                << "/preconditioner/near_null_space";// set a (or multiple) near null space(s)
  if (Spud::have_option(buffer.str()))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
    Mat pmat;
    perr = PCGetOperators(pc, PETSC_NULL, &pmat, PETSC_NULL);
    CHKERRV(perr);

    MatNullSpace SP;                                                 // create a set of nullspaces in a null space object
    fill_nullspace_(buffer.str(), SP, parent_indices);

    perr = MatSetNearNullSpace(pmat, SP); CHKERRV(perr);

    perr = MatNullSpaceDestroy(&SP); CHKERRV(perr);                  // destroy the null space object, necessary?
    #else
    dolfin::error("Can only set near null spaces with petsc > 3.2.");
    #endif
  }

  if ((preconditioner=="ml")||(preconditioner=="gamg"))
  {
    perr = PCSetUp(pc); CHKERRV(perr);                               // need to call this to prevent seg fault on kspview
  }                                                                  // BUT it has to happen after the near null space is set

}

//*******************************************************************|************************************************************//
// fill a pc object from the options tree assuming its a fieldsplit pc
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_pc_fieldsplit_(const std::string &optionpath, 
                                           PC &pc, const std::string prefix,
                                           const std::vector<uint>* parent_indices)
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code


  std::vector< std::vector<uint> > child_indices;                    // a vector of vectors to collect the child indices (the
                                                                     // subsets of the parent_indices vector (if associated) that
                                                                     // will themselves become the parent_indices on the next
                                                                     // recursion)
  std::vector<uint> prev_indices;

  buffer.str(""); buffer << optionpath << "/fieldsplit";
  int nsplits = Spud::option_count(buffer.str());                    // how many fieldsplits exist for this pc
  for (uint i = 0; i < nsplits; i++)                                 // loop over them all
  {
    buffer.str(""); buffer << optionpath << 
                                        "/fieldsplit[" << i << "]";
    std::vector<uint> indices;
    IS is;
    if (i==0)
    {
      fill_is_by_field_(buffer.str(), is,                            // setup an IS for each fieldsplit
                        indices, parent_indices,
                        NULL);

      prev_indices = indices;                                        // should already be sorted so don't bother doing it again
    }
    else
    {
      fill_is_by_field_(buffer.str(), is,                            // setup an IS for each fieldsplit
                        indices, parent_indices,
                        &prev_indices);

      prev_indices.insert(prev_indices.end(), indices.begin(), indices.end());
      std::sort(prev_indices.begin(), prev_indices.end());           // sort the vector of prev_indices
    }

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    buffer.str(""); buffer << optionpath << 
                                        "/fieldsplit[" 
                                         << i << "]/name";
    std::string fsname;
    serr = Spud::get_option(buffer.str(), fsname);
    spud_err(buffer.str(), serr);
    
    buffer.str(""); buffer << prefix << fsname;
    perr = PCFieldSplitSetIS(pc, buffer.str().c_str(), is);          // set the fs using that IS
    CHKERRV(perr);
    perr = ISDestroy(&is); CHKERRV(perr);                            // destroy the IS, necessary?
    #else
    perr = PCFieldSplitSetIS(pc, is); CHKERRV(perr);                 // set the fs using that IS
    perr = ISDestroy(is); CHKERRV(perr);                             // destroy the IS, necessary?
    #endif

    child_indices.push_back(indices);                                // record the indices of the global vector that made it into
  }                                                                  // IS (and hence the fieldsplit)

  std::string ctype;
  buffer.str(""); buffer << optionpath << "/composite_type/name";    // composite type of fieldsplit (additive, multiplicative etc.)
  serr = Spud::get_option(buffer.str(), ctype);                      // sadly no string based interface provided to this (I think)
  spud_err(buffer.str(), serr);                                      // so hard code an if block

  if (ctype == "additive")                                           // additive fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "multiplicative")                                // multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "symmetric_multiplicative")                      // symmetric multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, 
                          PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "special")                                       // special fieldsplit (whatever that means!)
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SPECIAL); 
    CHKERRV(perr);
  }
  else if (ctype == "schur")                                         // schur fieldsplit
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR); 
    CHKERRV(perr);
    
    std::string ftype;
    buffer.str(""); buffer << optionpath << 
                    "/composite_type::schur/factorization_type/name";// schur factorization type
    serr = Spud::get_option(buffer.str(), ftype);                    // sadly no string based interface provided to this (I think)
    spud_err(buffer.str(), serr);                                    // so hard code an if block
    if (ftype == "full")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
      perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL); 
      CHKERRV(perr);
      #endif
    }
    else
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
      if (ftype == "upper")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER); 
        CHKERRV(perr);
      }
      else if (ftype == "lower")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER); 
        CHKERRV(perr);
      }
      else if (ftype == "diag")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG); 
        CHKERRV(perr);
      }
      else
      {
        dolfin::error("Unknown PCFieldSplitSchurFactType.");
      }
      #else
      dolfin::error("Can only set schur factorization_type to anything other than full with petsc > 3.2.");
      #endif
    }
    
    std::string ptype;
    buffer.str(""); buffer << optionpath << 
                  "/composite_type::schur/schur_preconditioner/name";// preconditioner for schur block
    serr = Spud::get_option(buffer.str(), ptype);                    // sadly no string based interface provided to this (I think)
    spud_err(buffer.str(), serr);                                    // so hard code an if block
    if (ptype == "diag")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_A11, PETSC_NULL);
      #else
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_DIAG, PETSC_NULL);
      #endif
      CHKERRV(perr);
    }
    else if (ptype == "self")
    {
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_SELF, PETSC_NULL);
      CHKERRV(perr);
    }
    else if (ptype == "user")
    {
      PetscInt n = child_indices[1].size();
      PetscInt *petscindices;
      PetscMalloc(n*sizeof(PetscInt), &petscindices);
      
      uint ind = 0;
      for (std::vector<uint>::const_iterator                         // loop over the indices
                                  ind_it = child_indices[1].begin(); 
                                  ind_it != child_indices[1].end(); 
                                  ind_it++)
      {
        petscindices[ind] = *ind_it;                                 // insert them into the PetscInt array
        ind++;                                                       // increment the array index
      }
      assert(ind==n);                                                // these should be equal

      IS_ptr is;
      is.reset( new IS );
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      perr = ISCreateGeneral(PETSC_COMM_WORLD, n, petscindices, 
                                        PETSC_OWN_POINTER, &(*is));  // create the general index set based on the indices
      #else
      perr = ISCreateGeneral(PETSC_COMM_WORLD, n, 
                                               petscindices, &(*is));// create the general index set based on the indices
      #endif
      CHKERRV(perr);
      
      buffer.str(""); buffer << optionpath << 
                 "/composite_type::schur/schur_preconditioner::user";
      if (Spud::have_option(buffer.str()+"/monitors/view_index_set"))
      {
        std::string isname = prefix+"SchurPC";
        dolfin::log(dolfin::INFO, "ISView: %s (%s)", 
                             isname.c_str(), buffer.str().c_str());
        perr = ISView(*is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr); // isview?
      }

      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #else
      PetscFree(petscindices);                                       // free the PetscInt array of indices
      #endif

      solverindexsets_[prefix+"SchurPC"] = is;

      Mat_ptr submatrix;
      submatrix.reset( new Mat );
      perr = MatGetSubMatrix(*(*solvermatrices_[prefix+"SchurPC"]).mat(), *is, *is, MAT_INITIAL_MATRIX, &(*submatrix));
      CHKERRV(perr);

      solversubmatrices_[prefix+"SchurPC"] = submatrix;
      
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, *submatrix);
      CHKERRV(perr);
    }
    else
    {
      dolfin::error("Unknown PCFieldSplitSchurPreconditionType.");
    }
    #else
    dolfin::error("schur fieldsplits only support with petsc > 3.1");
    #endif
  }
  else                                                               // unknown (to buckettools) fieldsplit composite type
  {
    dolfin::error("Unknown PCCompositeType.");
  }
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  perr = PCSetUp(pc); CHKERRV(perr);                                 // call this before subksp can be retrieved
                                                                     // this is only necessary for schur fieldsplits and
                                                                     // as we're not supporting these with petsc 3.1 we don't
                                                                     // bother otherwise
  #endif

  KSP *subksps;                                                      // setup the fieldsplit subksps
  PetscInt nsubksps;
  perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); 
  CHKERRV(perr); 

  assert(nsubksps==nsplits);

  for (uint i = 0; i < nsplits; i++)                                 // loop over the splits again
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = KSPSetNormType(subksps[i], KSP_NORM_DEFAULT);             // these two calls are necessary because the normtype
    CHKERRV(perr);                                                   // and pc_side will have been set inappropriately in the 
    perr = KSPSetPCSide(subksps[i], PC_SIDE_DEFAULT);                // above call to PCSetUp.  now we undo this damage by resetting
    CHKERRV(perr);                                                   // them to default (since the schema doesn't allow them to be
                                                                     // set by the user anyway) so that they can once again be reset
                                                                     // by KSP/PCSetUp later (probably from KSPSolve) once the ksp
                                                                     // type is set properly.
                                                                     // same petsc 3.1 support logic as for PCSetUp above.
    #endif

    std::string fsname;
    buffer.str(""); buffer << optionpath << "/fieldsplit[" 
                                        << i << "]/name";
    serr = Spud::get_option(buffer.str(), fsname);
    spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << "/fieldsplit[" 
                                        << i << "]/linear_solver";
    fill_ksp_(buffer.str(), subksps[i], prefix+fsname+"_", 
                                                &child_indices[i]);  // recurse and fill in the data on each subksp but passing down
                                                                     // the child_indices as the new parent_indices
  }

}

//*******************************************************************|************************************************************//
// Fill a petsc IS object from the options tree (mostly for fieldsplits).
// IS's may be set up by field name, components of the field, regions of the domain of the field  and surfaces of the domain of the
// field.
// Additionally the resulting IS is checked for consistency with any parents or siblings in the tree.
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_is_by_field_(const std::string &optionpath, IS &is, 
                                         std::vector<uint> &child_indices, 
                                         const std::vector<uint>* parent_indices,
                                         const std::vector<uint>* sibling_indices)
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  buffer.str(""); buffer << optionpath << "/field";                  // loop over the fields used to describe this IS
  int nfields = Spud::option_count(buffer.str());

  child_indices.clear();
  if (nfields==0)                                                    // if no fields have been specified... **no fields**
  {
    std::pair<uint, uint> ownership_range =                          // the parallel ownership range of the system functionspace
            (*(*(*system_).functionspace()).dofmap()).ownership_range();
    for (uint i = ownership_range.first; i < ownership_range.second; i++)
    {
      child_indices.push_back(i);
    }
  }
  else
  {

    bool mixedsystem = (((*system_).fields_size())>1);

    for (uint i = 0; i < nfields; i++)                               // loop over the fields that have been specified
    {

      std::string fieldname;
      buffer.str(""); buffer << optionpath << "/field[" << i <<      // get the field name
                                                            "]/name";
      serr = Spud::get_option(buffer.str(), fieldname); 
      spud_err(buffer.str(), serr);
      
      FunctionBucket_ptr field = (*system_).fetch_field(fieldname);  // using the name, get the
      const int fieldindex = (*field).index();                       // field index
      const std::string fieldrank = (*field).rank();                 // field rank
      const int fieldsize = (*field).size();                         // field size

      FunctionSpace_ptr functionspace;
      if (mixedsystem)
      {
        functionspace = (*(*system_).functionspace())[fieldindex];
      }
      else
      {
        functionspace = (*system_).functionspace();
      }

      std::vector<int>* components = NULL;
      std::vector<int>* region_ids = NULL;
      std::vector<int>* boundary_ids = NULL;
      buffer.str(""); buffer << optionpath << "/field[" << i << "]";
      is_by_field_restrictions_(buffer.str(),                        // get any restrictions on the IS
                                components, region_ids, boundary_ids, 
                                fieldrank, fieldsize);
     
      boost::unordered_set<uint> dof_set;
      dof_set  = field_dof_set_(buffer.str(), functionspace, 
                                components, region_ids, boundary_ids);

      child_indices.insert(child_indices.end(), dof_set.begin(), dof_set.end());

      destroy_is_field_restrictions_(components, region_ids, boundary_ids);
    }
  }

  restrict_is_indices_(child_indices, parent_indices, 
                                      sibling_indices);

  PetscInt n=child_indices.size();                                   // setup a simpler structure for petsc
  assert(n>0);
  PetscInt *indices;
  PetscMalloc(n*sizeof(PetscInt), &indices);
 
  uint ind = 0;
  if(parent_indices)
  {                                                                  // we have been passed a list of parent indices... 
                                                                     // our child indices must be a  subset of this list and indexed
                                                                     // into it so let's do that now while we convert structures...
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        c_it = child_indices.begin(); 
                                        c_it != child_indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // child_indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and throw an error
          dolfin::error("IS indices are not a subset of a parent fieldsplit, shouldn't happen here.");
        }
      }
      indices[ind] = p_ind;                                          // found the child index in the parent_indices so copy it into
                                                                     // the PetscInt array
      ind++;                                                         // increment the array index
      p_ind++;                                                       // indices shouldn't be repeated so increment the parent too
    } 
    assert(ind==n);                                                  // these should be equal
  }
  else
  {
    for (std::vector<uint>::const_iterator                           // loop over the child_indices
                                      ind_it = child_indices.begin(); 
                                      ind_it != child_indices.end(); 
                                      ind_it++)
    {
      indices[ind] = *ind_it;                                        // insert them into the PetscInt array
      ind++;                                                         // increment the array index
    }
    // these should be equal
    assert(ind==n);                                                  // these should be equal
  }

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, 
                                    PETSC_OWN_POINTER, &is);         // create the general index set based on the indices
  #else
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, &is);         // create the general index set based on the indices
  #endif
  CHKERRV(perr);
  if (Spud::have_option(optionpath+"/monitors/view_index_set"))
  {
    buffer.str(""); buffer << optionpath << "/name";                 // IS Name
    std::string isname;
    serr = Spud::get_option(buffer.str(), isname);
    spud_err(buffer.str(), serr);
    
    dolfin::log(dolfin::INFO, "ISView: %s (%s)", 
                                isname.c_str(), optionpath.c_str());
    perr = ISView(is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);      // isview?
  }

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  #else
  PetscFree(indices);                                                // free the PetscInt array of indices
  #endif
   
}

//*******************************************************************|************************************************************//
// return a vector of dofs from the given functionspace for a field
//*******************************************************************|************************************************************//
boost::unordered_set<uint> SpudSolverBucket::field_dof_set_(const std::string &optionpath,
                                                            const FunctionSpace_ptr functionspace,
                                                            const std::vector<int>* components,
                                                            const std::vector<int>* region_ids,
                                                            const std::vector<int>* boundary_ids,
                                                            const uint parent_component,
                                                            uint rank)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  assert(rank<=2);                                                   // component logic below only makes sense for rank <= 2

  boost::unordered_set<uint> dof_set;

  const uint num_sub_elements = (*(*functionspace).element()).num_sub_elements();
  if (num_sub_elements>0)
  {

    for (uint i = 0; i < num_sub_elements; i++)
    {
      if (components)
      {
        std::vector<int>::const_iterator comp = std::find((*components).begin(), 
                                           (*components).end(), 
                                           parent_component*num_sub_elements + i);
        if (comp == (*components).end())
        {
          continue;                                                  // component not requested so continue
        }
      }

      boost::unordered_set<uint> tmp_dof_set;
      tmp_dof_set = field_dof_set_(optionpath, (*functionspace)[i], 
                                   components, region_ids, boundary_ids,
                                   (parent_component*num_sub_elements + i), 
                                   rank++);
      dof_set.insert(tmp_dof_set.begin(), tmp_dof_set.end());
    }

    return dof_set;
  }
  
  assert(num_sub_elements==0);

  boost::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();

  if (boundary_ids)                                                  // do we have boundary id restrictions
  {                                                                  // yes, then get the dofs over these boundaries
    if (region_ids)
    {
      dof_set = cell_dof_set_(dofmap, region_ids);                   // if we have boundary_ids then we're only interested
    }                                                                // in cell dofs if we have region_ids specified too
    boost::unordered_set<uint> facet_dof_set;
    facet_dof_set = facet_dof_set_(dofmap, boundary_ids);
    dof_set.insert(facet_dof_set.begin(), facet_dof_set.end());
  }
  else                                                               // no boundary_ids specified so let's hope we have some
  {                                                                  // cells to fill the goody bag with
    dof_set = cell_dof_set_(dofmap, region_ids);
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// return a vector of dofs from the given dofmap possibly for a subset of the region ids as specified
//*******************************************************************|************************************************************//
boost::unordered_set<uint> SpudSolverBucket::cell_dof_set_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                           const std::vector<int>* region_ids)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  boost::unordered_set<uint> dof_set;

  Mesh_ptr mesh = (*system_).mesh();                                 // get the mesh
  MeshFunction_size_t_ptr cellidmeshfunction;

  if (region_ids)
  {                                                                  // yes...  **field(+component)+region(+boundary)**
    cellidmeshfunction = (*system_).celldomains();
  }

  for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)        // loop over the cells in the mesh
  {
    if (region_ids)
    {
      int cellid = (*cellidmeshfunction)[(*cell).index()];           // get the cell region id from the mesh function

      std::vector<int>::const_iterator id = std::find((*region_ids).begin(), 
                                             (*region_ids).end(), cellid);
      if (id == (*region_ids).end())
      {
        continue;
      }
    }

    std::vector<dolfin::la_index> dof_vec = (*dofmap).cell_dofs((*cell).index());
    for (std::vector<dolfin::la_index>::const_iterator dof_it =           // loop over the cell dof
                                    dof_vec.begin(); 
                                    dof_it < dof_vec.end(); 
                                    dof_it++)
    {
      dof_set.insert((uint) *dof_it);                                       // and insert each one into the unordered set
    }                                                                // (i.e. if it hasn't been added already)
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// return a vector of dofs from the given dofmap for the boundary ids specified
//*******************************************************************|************************************************************//
boost::unordered_set<uint> SpudSolverBucket::facet_dof_set_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                            const std::vector<int>* boundary_ids)
{
  Spud::OptionError serr;                                            // spud error code

  boost::unordered_set<uint> dof_set;                                // set up an unordered set of dof

  Mesh_ptr mesh = (*system_).mesh();                                 // get the mesh
  MeshFunction_size_t_ptr facetidmeshfunction = (*system_).facetdomains();

  for (dolfin::FacetIterator facet(*mesh); !facet.end(); ++facet)    // loop over the facets in the mesh
  {
    int facetid = (*facetidmeshfunction)[(*facet).index()];          // get the facet region id from the mesh function

    for (std::vector<int>::const_iterator id =                       // loop over the region ids that have been requested
                                (*boundary_ids).begin(); 
                                id != (*boundary_ids).end(); id++)
    {
      if(facetid==*id)                                               // check if this facet should be included
      {                                                              // yes...

        const dolfin::Cell cell(*mesh,                               // get cell to which facet belongs
               (*facet).entities((*mesh).topology().dim())[0]);      // (there may be two, but pick first)

        const std::size_t facet_number = cell.index(*facet);                // get the local index of the facet w.r.t. the cell

        std::vector<dolfin::la_index> cell_dof_vec;
        cell_dof_vec = (*dofmap).cell_dofs(cell.index());            // get the cell dof (potentially for all components)
        
        std::vector<std::size_t> facet_dof_vec((*dofmap).num_facet_dofs(), 0);
        (*dofmap).tabulate_facet_dofs(facet_dof_vec, facet_number);

        for (std::vector<std::size_t>::const_iterator dof_it =              // loop over the cell dof
                                facet_dof_vec.begin(); 
                                dof_it < facet_dof_vec.end(); 
                                dof_it++)
        {
          dof_set.insert((uint) cell_dof_vec[*dof_it]);                     // and insert each one into the unordered set
        }                                                            // (i.e. if it hasn't been added already)
      }
    }
  }

  return dof_set;

}

//*******************************************************************|************************************************************//
// fill a petsc nullspace object using the options in the optionpath provided
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_nullspace_(const std::string &optionpath, MatNullSpace &SP,
                                       const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  buffer.str(""); buffer << optionpath << "/null_space";
  int nnulls = Spud::option_count(buffer.str());                     // how many null spaces?

  std::vector< PETScVector_ptr > nullvecs;                           // collect the null space vectors here (so we maintain a reference)
  Vec vecs[nnulls];                                                  // and here (for the petsc interface)

  uint kspsize = 0;
  if(parent_indices)
  {                                                                  // if parent_indices is associated then this is a null space
    kspsize = (*parent_indices).size();                              // of a subksp so the kspsize is not the whole thing
  }                                                                  // FIXME: broken in parallel?
  else
  {                                                                 
    kspsize = (*(*(*system_).function()).vector()).size();           // otherwise, this is quite easy - just the size of the parent
  }                                                                  // system function

  for (uint i = 0; i<nnulls; i++)                                    // loop over the nullspaces
  {

    PETScVector_ptr nullvec( new dolfin::PETScVector(kspsize, "local") );// create a null vector for this null space

    buffer.str(""); buffer << optionpath <<                          // optionpath of the nullspace
                      "/null_space[" << i << "]";
    fill_values_by_field_(buffer.str(), nullvec, 0.0,                // create a vector describing the nullspace based on this optionpath 
                                parent_indices, NULL);               // (no siblings as null spaces can overlap)

    perr = VecNormalize(*(*nullvec).vec(), PETSC_NULL); CHKERRV(perr);// normalize the null space vector

    nullvecs.push_back(nullvec);                                     // keep the null vector in scope by grabbing a reference to it
    vecs[i] = *(*nullvec).vec();                                     // also collect it in a petsc compatible format (shouldn't take
                                                                     // reference though... hence line above, necessary?)
  
  }

  perr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, nnulls, 
                                                      vecs, &SP); 
  CHKERRV(perr);

  buffer.str(""); buffer << optionpath << 
                    "/monitors/view_null_space";                     // view the null space for debugging
  if (Spud::have_option(buffer.str()))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatNullSpaceView(SP, PETSC_VIEWER_STDOUT_SELF); 
    CHKERRV(perr);
    #else
    dolfin::log(dolfin::WARNING, "Cannot set view_null_space monitor with PETSc < 3.2.");
    #endif
  }
}

//*******************************************************************|************************************************************//
// fill any constraints used on fields using the optionpath provided
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_constraints_()
{
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  std::stringstream buffer;                                          // optionpath buffer
  PetscErrorCode perr;                                               // petsc error code

  PETScVector_ptr ub;
  buffer.str(""); buffer << optionpath() << "/type/snes_type/constraints/upper_bound";
  fill_bound_(buffer.str(), ub, SNES_VI_INF);

  PETScVector_ptr lb;
  buffer.str(""); buffer << optionpath() << "/type/snes_type/constraints/lower_bound";
  fill_bound_(buffer.str(), lb, SNES_VI_NINF);

  perr = SNESVISetVariableBounds(snes_, *(*lb).vec(), *(*ub).vec());
  CHKERRV(perr);
                                                                     // UGLY HACK: our constant bounds will be overwritten by the dm
                                                                     // in SNESSetUp so let's stop it from doing that by attaching a
                                                                     // dummy variable bounds computation - does nothing!
  perr = SNESVISetComputeVariableBounds(snes_, SNESVIDummyComputeVariableBounds);
  CHKERRV(perr);
  #else
  dolfin::error("SNES VI only available with petsc > 3.1");
  #endif

}

//*******************************************************************|************************************************************//
// fill petsc vectors describing the bounds on fields using the optionpath provided
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_bound_(const std::string &optionpath, PETScVector_ptr &bound, const double &background_value)
{

  uint size = 0;
  size = (*(*(*system_).function()).vector()).local_size();
  bound.reset( new dolfin::PETScVector(size, "local") );

  if (Spud::have_option(optionpath))
  {
    fill_values_by_field_(optionpath, bound, background_value,        // create a vector describing a bound based on this optionpath 
                                NULL, NULL);
  }
  else
  {
    std::vector<double> background(size, background_value);
    (*bound).set_local(background);
    background.clear();
  }

}

//*******************************************************************|************************************************************//
// Fill a vector of values from the options tree.
// IS's may be set up by field name, components of the field, regions of the domain of the field  and surfaces of the domain of the
// field.
// Additionally the resulting IS is checked for consistency with any parents in the tree.
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_values_by_field_(const std::string &optionpath, PETScVector_ptr values,
                                             const double &background_value,
                                             const std::vector<uint>* parent_indices,
                                             const std::vector<uint>* sibling_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  buffer.str(""); buffer << optionpath << "/field";                  // loop over the fields used to describe this vector
  int nfields = Spud::option_count(buffer.str());

  boost::unordered_map<uint, double> value_map;
  if (nfields==0)                                                    // if no fields have been specified...
  {
    std::pair<uint, uint> ownership_range =                          // the parallel ownership range of the system functionspace
            (*(*(*system_).functionspace()).dofmap()).ownership_range();
    for (uint i = ownership_range.first; i < ownership_range.second; i++)
    {
      value_map[i] = 1.0;                                            // we assume a constant value of one
    }
  }
  else
  {

    bool mixedsystem = (((*system_).fields_size())>1);

    for (uint i = 0; i < nfields; i++)                               // loop over the fields that have been specified
    {

      std::string fieldname;
      buffer.str(""); buffer << optionpath << "/field[" << i <<      // get the field name
                                                            "]/name";
      serr = Spud::get_option(buffer.str(), fieldname); 
      spud_err(buffer.str(), serr);
      
      FunctionBucket_ptr field = (*system_).fetch_field(fieldname);  // using the name, get the
      const int fieldindex = (*field).index();                       // field index
      const std::string fieldrank = (*field).rank();                 // field rank
      const int fieldsize = (*field).size();                         // and field size

      FunctionSpace_ptr functionspace;                               // grab the functionspace
      if (mixedsystem)
      {
        functionspace = (*(*system_).functionspace())[fieldindex];
      }
      else
      {
        functionspace = (*system_).functionspace();
      }

      std::vector<int>* components = NULL;
      std::vector<int>* region_ids = NULL;
      std::vector<int>* boundary_ids = NULL;
      buffer.str(""); buffer << optionpath 
                                    << "/field[" << i << "]";
      is_by_field_restrictions_(buffer.str(),                        // get the standard restrictions on an IS
                                components, region_ids, boundary_ids,
                                fieldrank, fieldsize);
     
      Expression_ptr value_exp;
      buffer.str(""); buffer << optionpath << "/field[" 
                                            << i << "]/python";
      if (Spud::have_option(buffer.str()))
      {
        std::string pyfunction;
        serr = Spud::get_option(buffer.str(), pyfunction);
        spud_err(buffer.str(), serr);

        if (fieldrank=="Scalar")
        {
          value_exp.reset( new PythonExpression(pyfunction) );
        }
        else if (fieldrank=="Vector")
        {
          int size = fieldsize;
          if (components)
          {
            size = (*components).size();
          }
          value_exp.reset( new PythonExpression(size, pyfunction) );
        }
        else
        {
          dolfin::error("Tensor value maps not implemented yet.");
        }
      }

      double *value_const = NULL;
      buffer.str(""); buffer << optionpath << "/field[" 
                                            << i << "]/constant";
      if (Spud::have_option(buffer.str()))
      {
        value_const = new double;
        serr = Spud::get_option(buffer.str(), *value_const);
        spud_err(buffer.str(), serr);
      }

      boost::unordered_map<uint, double> tmp_map;
      buffer.str(""); buffer << optionpath << "/field[" << i << "]";
      tmp_map = field_value_map_(buffer.str(), functionspace,        // for each field construct a map describing the indices and
                             components, region_ids,                 // values
                             boundary_ids, value_exp, value_const);
      value_map.insert(tmp_map.begin(), tmp_map.end());

      destroy_is_field_restrictions_(components, 
                                     region_ids, 
                                     boundary_ids);
      if(value_const)
      {
        delete(value_const);
        value_const = NULL;
      }

    }
  }

  std::vector<uint> child_indices;                    
  for (boost::unordered_map<uint, double>::const_iterator            // construct a vector of indices which we'll restrict and order then use
                                              c_it = value_map.begin(); // to index the map
                                              c_it != value_map.end();
                                              c_it++)
  {
    child_indices.push_back((*c_it).first);
  }

  restrict_is_indices_(child_indices, parent_indices, sibling_indices);// restrict based on the parents (intersection) and possibly
                                                                     // siblings 

  PetscInt n=child_indices.size();                                   // setup a simpler structure for petsc
  assert(n>0);
  PetscInt *indices;
  PetscMalloc(n*sizeof(PetscInt), &indices);
  dolfin::PETScVector vvec(n, "local");                             // create a local vector of local size length 
 
  dolfin::la_index ind = 0;
  if(parent_indices)
  {                                                                  // we have been passed a list of parent indices... 
                                                                     // our child indices must be a  subset of this list and indexed
                                                                     // into it so let's do that now while we convert structures...
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        c_it = child_indices.begin(); 
                                        c_it != child_indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // child_indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and throw an error
          dolfin::error("IS indices are not a subset of a parent fieldsplit, shouldn't happen here.");
        }
      }
      indices[ind] = p_ind;                                          // found the child index in the parent_indices so copy it into
                                                                     // the PetscInt array
      vvec.set(&value_map[(*c_it)], 1, &ind);                          // set the null vector to the value in the map
      ind++;                                                         // increment the array index
      p_ind++;                                                       // indices shouldn't be repeated so increment the parent too
    } 
    assert(ind==n);                                                  // these should be equal
  }
  else
  {
    for (std::vector<uint>::const_iterator                           // loop over the child_indices
                                      ind_it = child_indices.begin(); 
                                      ind_it != child_indices.end(); 
                                      ind_it++)
    {
      indices[ind] = *ind_it;                                        // insert them into the PetscInt array
      vvec.set(&value_map[(*ind_it)], 1, &ind);                        // set the null vector to the value in the map
      ind++;                                                         // increment the array index
    }
    // these should be equal
    assert(ind==n);                                                  // these should be equal
  }

  IS is;
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, 
                                    PETSC_OWN_POINTER, &is);         // create the general index set based on the indices
  #else
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, &is);         // create the general index set based on the indices
  PetscFree(indices);                                                // free the PetscInt array of indices
  #endif
  CHKERRV(perr);
   
  if (Spud::have_option(optionpath+"/monitors/view_index_set"))
  {
    buffer.str(""); buffer << optionpath << "/name";                 // IS Name
    if (Spud::have_option(buffer.str()))
    {
      std::string isname;
      serr = Spud::get_option(buffer.str(), isname);
      spud_err(buffer.str(), serr);
      dolfin::log(dolfin::INFO, "ISView: %s (%s)", 
                                  isname.c_str(), optionpath.c_str());
    }
    else
    {
      dolfin::log(dolfin::INFO, "ISView: (%s)", 
                                  optionpath.c_str());
    }
    perr = ISView(is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);      // isview?
  }

  std::vector<double> background((*values).local_size(), background_value);
  (*values).set_local(background);
  background.clear();

  VecScatter scatter;                                                // create a petsc scatter object from an object with the same 
  perr = VecScatterCreate(*vvec.vec(), PETSC_NULL,                  // structure as the vvec vector to one with the same structure
                          *(*values).vec(), is, &scatter);          // as the null vector using the IS
  CHKERRV(perr);
  perr = VecScatterBegin(scatter, *vvec.vec(),                      // scatter from the vvec vector to the null vector
                         *(*values).vec(), INSERT_VALUES, 
                         SCATTER_FORWARD); 
  CHKERRV(perr);
  perr = VecScatterEnd(scatter, *vvec.vec(), 
                       *(*values).vec(), INSERT_VALUES, 
                       SCATTER_FORWARD); 
  CHKERRV(perr);
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1            // necessary or taken care of when object leaves scope?
  perr = VecScatterDestroy(&scatter); CHKERRV(perr);      
  perr = ISDestroy(&is); CHKERRV(perr);
  #else
  perr = VecScatterDestroy(scatter); CHKERRV(perr);
  perr = ISDestroy(is); CHKERRV(perr);
  #endif

  (*values).apply("");                                              // finish assembly of the null vector, just in case

}

//*******************************************************************|************************************************************//
// return a map from dofs to null space values from the given functionspace for a field
//*******************************************************************|************************************************************//
boost::unordered_map<uint, double> SpudSolverBucket::field_value_map_(const std::string &optionpath,
                                                                   const FunctionSpace_ptr functionspace,
                                                                   const std::vector<int>* components,
                                                                   const std::vector<int>* region_ids,
                                                                   const std::vector<int>* boundary_ids, 
                                                                   Expression_ptr value_exp, const double *value_const,
                                                                   const uint parent_component,
                                                                   uint rank, uint exp_index)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  assert(rank<=2);                                                   // component logic below only makes sense for rank <= 2

  boost::unordered_map<uint, double> value_map;

  const uint num_sub_elements = (*(*functionspace).element()).num_sub_elements();

  if (num_sub_elements>0)
  {

    for (uint i = 0; i < num_sub_elements; i++)
    {
      if (components)
      {
        std::vector<int>::const_iterator comp = std::find((*components).begin(), 
                                           (*components).end(), 
                                           parent_component*num_sub_elements + i);
        if (comp == (*components).end())
        {
          continue;                                                  // component not requested so continue
        }
        exp_index = comp - (*components).begin();                    // work out the index into the expression for this component
      }
      else
      {
        exp_index = i;
      }

      boost::unordered_map<uint, double> tmp_map;
      tmp_map = field_value_map_(optionpath, (*functionspace)[i], 
                                 components, region_ids,
                                 boundary_ids, value_exp, value_const, 
                                 (parent_component*num_sub_elements + i), 
                                 rank++, exp_index);
      value_map.insert(tmp_map.begin(), tmp_map.end());
    }

    return value_map;
  }
  
  assert(num_sub_elements==0);

  boost::shared_ptr<const dolfin::GenericDofMap> dofmap = (*functionspace).dofmap();

  if (boundary_ids)                                                  // do we have boundary ids specified?
  {                                                                  // yes...
    if (region_ids)                                                  // if we have boundary_ids then we're only interested
    {                                                                // in cell dofs if we have region_ids specified too
      value_map = cell_value_map_(dofmap, region_ids, value_exp, value_const, exp_index);
    }
    boost::unordered_map<uint, double> facet_value_map;
    facet_value_map = facet_value_map_(dofmap, boundary_ids, value_exp, value_const, exp_index);
    value_map.insert(facet_value_map.begin(), facet_value_map.end());
  }
  else
  {
    value_map = cell_value_map_(dofmap, region_ids, value_exp, value_const, exp_index);
  }

  return value_map;

}

//*******************************************************************|************************************************************//
// return an unordered map from dofs to null space values from the given dofmap possibly just for the region ids specified
//*******************************************************************|************************************************************//
boost::unordered_map<uint, double> SpudSolverBucket::cell_value_map_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                                  const std::vector<int>* region_ids,
                                                                  Expression_ptr value_exp, const double* value_const,
                                                                  const uint &exp_index)
{
  Spud::OptionError serr;                                            // spud error code

  boost::unordered_map<uint, double> value_map;

  Mesh_ptr mesh = (*system_).mesh();                                 // get the mesh

  const uint gdim = (*mesh).geometry().dim();                        // set up data for expression evaluation
  boost::multi_array<double, 2> coordinates(boost::extents[(*dofmap).max_cell_dimension()][gdim]);
  std::vector<double> vertex_coordinates;
  dolfin::Array<double> x(gdim);
  uint value_size = 1;
  if (value_exp)
  {
    for (uint i = 0; i < (*value_exp).value_rank(); i++)
    {
      value_size *= (*value_exp).value_dimension(i);
    }
    assert(!value_const);
  }
  else
  {
    assert(value_const);
  }
  dolfin::Array<double> values(value_size);

  MeshFunction_size_t_ptr cellidmeshfunction;
  if (region_ids)
  {
    cellidmeshfunction = (*system_).celldomains();
  }

  for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)        // loop over the cells in the mesh
  {
    if (region_ids)
    {
      int cellid = (*cellidmeshfunction)[(*cell).index()];           // get the cell region id from the mesh function

      std::vector<int>::const_iterator id = std::find((*region_ids).begin(), 
                                             (*region_ids).end(), cellid);
      if (id == (*region_ids).end())
      {
        continue;                                                    // region id not requested so continue
      }
    }

    std::vector<dolfin::la_index> dof_vec = (*dofmap).cell_dofs((*cell).index());

    if(value_exp)
    {
      (*cell).get_vertex_coordinates(vertex_coordinates);
      (*dofmap).tabulate_coordinates(coordinates, vertex_coordinates, *cell);
    }

    for (uint i = 0; i < dof_vec.size(); i++)                        // loop over the cell dof
    {
      if(value_exp)
      {
        for (uint j = 0; j < gdim; j++)
        {
          x[j] = coordinates[i][j];
        }
        (*value_exp).eval(values, x);                                // evaluate te expression
        value_map[(uint) dof_vec[i]] = values[exp_index];            // and set the null space to that
      }
      else
      {
        value_map[(uint) dof_vec[i]] = *value_const;                 // and insert each one into the unordered map
                                                                     // assuming a constant
      }
    }
  }

  return value_map;

}

//*******************************************************************|************************************************************//
// return an unordered map from dofs to null space values from the given dofmap for the boundary ids specified
//*******************************************************************|************************************************************//
boost::unordered_map<uint, double> SpudSolverBucket::facet_value_map_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                                   const std::vector<int>* boundary_ids, 
                                                                   Expression_ptr value_exp, const double *value_const, 
                                                                   const uint &exp_index)
{
  Spud::OptionError serr;                                            // spud error code

  assert(boundary_ids);

  boost::unordered_map<uint, double> value_map;                         // set up an unordered set of dof

  Mesh_ptr mesh = (*system_).mesh();                                 // get the mesh

  const uint gdim = (*mesh).geometry().dim();
  boost::multi_array<double, 2> coordinates(boost::extents[(*dofmap).max_cell_dimension()][gdim]);
  std::vector<double> vertex_coordinates;
  dolfin::Array<double> x(gdim);
  uint value_size = 1;
  if (value_exp)
  {
    for (uint i = 0; i < (*value_exp).value_rank(); i++)
    {
      value_size *= (*value_exp).value_dimension(i);
    }
    assert(!value_const);
  }
  else
  {
    assert(value_const);
  }
  dolfin::Array<double> values(value_size);

  MeshFunction_size_t_ptr facetidmeshfunction = (*system_).facetdomains();

  for (dolfin::FacetIterator facet(*mesh); !facet.end(); ++facet)    // loop over the facets in the mesh
  {
    int facetid = (*facetidmeshfunction)[(*facet).index()];          // get the facet boundary id from the mesh function

    for (std::vector<int>::const_iterator id =                       // loop over the boundary ids that have been requested
                                (*boundary_ids).begin(); 
                                id != (*boundary_ids).end(); id++)
    {
      if(facetid==*id)                                               // check if this facet should be included
      {                                                              // yes...

        const dolfin::Cell cell(*mesh,                               // get cell to which facet belongs
               (*facet).entities((*mesh).topology().dim())[0]);      // (there may be two, but pick first)

        const std::size_t facet_number = cell.index(*facet);         // get the local index of the facet w.r.t. the cell

        std::vector<dolfin::la_index> cell_dof_vec;
        cell_dof_vec = (*dofmap).cell_dofs(cell.index());            // get the cell dof (potentially for all components)
        
        std::vector<std::size_t> facet_dof_vec((*dofmap).num_facet_dofs(), 0);
        (*dofmap).tabulate_facet_dofs(facet_dof_vec, facet_number);

        if (value_exp)
        {
          cell.get_vertex_coordinates(vertex_coordinates);
          (*dofmap).tabulate_coordinates(coordinates, vertex_coordinates, cell);
        }

        for (uint i = 0; i < facet_dof_vec.size(); i++)              // loop over facet dof
        {
          if(value_exp)
          {
            for (uint j = 0; j < gdim; j++)
            {
              x[j] = coordinates[i][j];
            }
            (*value_exp).eval(values, x);                            // evaluate the null space expression
            value_map[(uint) cell_dof_vec[facet_dof_vec[i]]] = values[exp_index];
          }
          else
          {
            value_map[(uint) cell_dof_vec[facet_dof_vec[i]]] = *value_const;// and insert each one into the unordered map
          }                                                          // assuming a constant
        }                                                         
      }
    }
  }

  return value_map;

}

//*******************************************************************|************************************************************//
// using the optionpath set up any restrictions we place on the field IS (i.e. region ids, components or boundary_ids)
//*******************************************************************|************************************************************//
void SpudSolverBucket::is_by_field_restrictions_(const std::string &optionpath,
                                                 std::vector<int>* &components,
                                                 std::vector<int>* &region_ids,
                                                 std::vector<int>* &boundary_ids,
                                                 const std::string &fieldrank,
                                                 const int &fieldsize)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath << "/components";
  if (Spud::have_option(buffer.str()))
  {                                            
    components = new std::vector<int>;
    serr = Spud::get_option(buffer.str(), *components);            // get the components
    spud_err(buffer.str(), serr);

    if(fieldrank=="Scalar")
    {
      dolfin::error("Requested null space components of a scalar field.");
    }
    else if(fieldrank=="Vector")
    {
      std::vector<int>::iterator max_comp_it =                  
           std::max_element((*components).begin(), (*components).end()); // check the maximum requested component exists
      
      assert(*max_comp_it < fieldsize);
      assert((*components).size() <= fieldsize);
    }
    else
    {
      dolfin::error("Only deal with scalar and vector null spaces for now.");
    }
  }

  buffer.str(""); buffer << optionpath << "/region_ids";
  if (Spud::have_option(buffer.str()))
  {                                             
    region_ids = new std::vector<int>;
    serr = Spud::get_option(buffer.str(), *region_ids);                // get the region ids
    spud_err(buffer.str(), serr);
  }

  buffer.str(""); buffer << optionpath << "/boundary_ids";
  if (Spud::have_option(buffer.str()))
  {                                              
    boundary_ids = new std::vector<int>;
    serr = Spud::get_option(buffer.str(), *boundary_ids);                // get the boundary ids
    spud_err(buffer.str(), serr);
  }


}

//*******************************************************************|************************************************************//
// destroy any restrictions we place on the field IS (i.e. region ids, components or boundary_ids)
//*******************************************************************|************************************************************//
void SpudSolverBucket::destroy_is_field_restrictions_(std::vector<int>* &components,
                                                      std::vector<int>* &region_ids,
                                                      std::vector<int>* &boundary_ids)
{
  if(components)
  {
    delete components;
    components = NULL;
  }
  if (region_ids)
  {
    delete region_ids;
    region_ids = NULL;
  }
  if (boundary_ids)
  {
    delete boundary_ids;
    boundary_ids = NULL;
  }

}

//*******************************************************************|************************************************************//
// restrict a vector of indices by its parent (intersection) or sibling (complement), also by parallel ownership
//*******************************************************************|************************************************************//
void SpudSolverBucket::restrict_is_indices_(std::vector<uint> &indices, 
                                            const std::vector<uint>* parent_indices, 
                                            const std::vector<uint>* sibling_indices)
{

  std::vector<uint> tmp_indices = indices;
  indices.clear();

  std::pair<uint, uint> ownership_range =                            // the parallel ownership range of the system functionspace
          (*(*(*system_).functionspace()).dofmap()).ownership_range();

  for (std::vector<uint>::const_iterator                             // loop over the dof in the set
                        dof_it = tmp_indices.begin(); 
                        dof_it != tmp_indices.end(); 
                        dof_it++)
  {                                                                  // and insert them into the indices vector
    if ((*dof_it >= ownership_range.first) &&                        // but first check that this process owns them
                          (*dof_it < ownership_range.second))        // (in parallel)
    {
      indices.push_back(*dof_it);
    }
  }

  std::sort(indices.begin(), indices.end());                         // sort the vector of indices

  if(sibling_indices)                                                // we have been passed a list of sibling indices...
  {                                                                  // we wish to remove from the indices any indices that
                                                                     // also occur in the sibling indices
    tmp_indices.clear();

    uint c_size = indices.size();
    uint c_ind = 0;
    bool overlap = false;
    for(std::vector<uint>::const_iterator                            // loop over the sibling indices
                                   s_it = (*sibling_indices).begin();
                                   s_it != (*sibling_indices).end();
                                   s_it++)
    {
      while(indices[c_ind] != *s_it)                                 // indices are sorted, so sibling_indices should be too
      {                                                              // search indices until the current sibling index is found
        tmp_indices.push_back(indices[c_ind]);                       // include indices that aren't in the sibling
        c_ind++;
        if (c_ind == c_size)                                         // or we reach the end of the indices...
        {
          break;
        }
      }
      if (c_ind == c_size)                                           // we've reached the end of the child indices so nothing more
      {                                                              // to do
        break;
      }
      else                                                           // we haven't reached the end of the child indices but found
      {                                                              // a sibling index to ignore... give a warning
        overlap = true;
      }
      c_ind++;                                                       // indices shouldn't be repeated so incredment the child too
    }
    for (uint i = c_ind; i < c_size; i++)
    {
      tmp_indices.push_back(indices[i]);                             // insert any remaining indices beyond the siblings
    }

    if(overlap)
    {                                                                // sibling indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices overlap with sibling fieldsplit, ignoring overlapping indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }

  if(parent_indices)                                                 // we have been passed a list of parent indices... 
  {                                                                  // we wish to remove from the indices any indices that do
                                                                     // not occur in the parent indices 
    tmp_indices.clear();

    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    uint p_reset = 0;
    bool extra = false;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        c_it = indices.begin(); 
                                        c_it != indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and prepare to throw a warning
          extra = true;
          break;
        }
      }
      if (p_ind == p_size)
      {
        p_ind = p_reset;
      }
      else
      {
        tmp_indices.push_back(*c_it);                                // include indices that are in the parent
        p_ind++;                                                     // indices shouldn't be repeated so increment the parent too
        p_reset = p_ind;                                             // this is where the next failed search should continue from
        if (p_ind == p_size)                                         // we've reached the end
        { 
          break;                                            
        }
      }
    } 

    if(extra)
    {                                                                // child indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices not a subset of parent fieldsplit, ignoring extra indices.");
    }
    indices.clear();
    indices = tmp_indices;
                                 
  }
}

//*******************************************************************|************************************************************//
// empty the data structures in the spudsolver bucket
//*******************************************************************|************************************************************//
void SpudSolverBucket::empty_()
{
  form_optionpaths_.clear();
  SolverBucket::empty_();
}

