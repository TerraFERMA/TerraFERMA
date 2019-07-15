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


#include "PythonExpression.h"
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
#include "DolfinPETScBase.h"
#include "Logger.h"
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

    perr = SNESCreate((*(*system_).mesh()).mpi_comm(), &snes_); 
    petsc_err(perr);                                                   // create the petsc snes object

    perr = SNESSetOptionsPrefix(snes_, prefix.str().c_str());        // set its petsc options name prefix to SystemName_SolverName
    petsc_err(perr);

    ctx_.solver = this;                                              // the snes context just needs this class... neat, huh?

    perr = SNESSetFunction(snes_, (*res_).vec(),                    // set the snes function to use the newly allocated residual vector
                                    FormFunction, (void *) &ctx_); 
    petsc_err(perr);

    if (bilinearpc_)                                                 // if we have a pc form
    {
      assert(matrixpc_);
      perr = SNESSetJacobian(snes_, (*matrix_).mat(),               // set the snes jacobian to have two matrices
                  (*matrixpc_).mat(), FormJacobian, (void *) &ctx_); 
      petsc_err(perr);
    }
    else                                                             // otherwise
    {
      perr = SNESSetJacobian(snes_, (*matrix_).mat(),               // set the snes jacobian to have the same matrix twice
                    (*matrix_).mat(), FormJacobian, (void *) &ctx_); 
      petsc_err(perr);
    }

    std::string snestype;
    buffer.str(""); buffer << optionpath() << "/type/snes_type/name";
    serr = Spud::get_option(buffer.str(), snestype);                 // set the snes type... ls is most common
    spud_err(buffer.str(), serr);

    if(snestype=="vi")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #if PETSC_VERSION_MINOR > 3
      perr = SNESSetType(snes_, SNESVINEWTONRSLS); petsc_err(perr); 
      #elif PETSC_VERSION_MINOR == 3
      perr = SNESSetType(snes_, SNESVIRS); petsc_err(perr); 
      #else
      perr = SNESSetType(snes_, snestype.c_str()); petsc_err(perr); 
      #endif
      perr = SNESSetFromOptions(snes_); petsc_err(perr);               // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)

      fill_constraints_();
      #else
      tf_err("Cannot set snes vi with PETSc < 3.2.", "PETSc version: %d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
      #endif
    }
    else if(snestype=="ls")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
      perr = SNESSetType(snes_, SNESNEWTONLS); petsc_err(perr);
      #else
      perr = SNESSetType(snes_, snestype.c_str()); petsc_err(perr);
      #endif 
      perr = SNESSetFromOptions(snes_); petsc_err(perr);               // set-up snes from options (we do this first to ensure that
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
        perr = SNESGetLineSearch(snes_, &linesearch); petsc_err(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); petsc_err(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "bt"); petsc_err(perr);
        perr = SNESLineSearchSetOrder(linesearch, 3); petsc_err(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchCubic, PETSC_NULL); petsc_err(perr); 
        #endif
      }
      else if (lstype=="quadratic")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); petsc_err(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); petsc_err(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "bt"); petsc_err(perr);
        perr = SNESLineSearchSetOrder(linesearch, 2); petsc_err(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchQuadratic, PETSC_NULL); petsc_err(perr); 
        #endif
      }
      else if (lstype=="basic")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); petsc_err(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); petsc_err(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "basic"); petsc_err(perr);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchNo, PETSC_NULL); petsc_err(perr); 
        #endif
      }
      else if (lstype=="basicnonorms")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
        SNESLineSearch linesearch;
        #if PETSC_VERSION_MINOR > 3
        perr = SNESGetLineSearch(snes_, &linesearch); petsc_err(perr);
        #else
        perr = SNESGetSNESLineSearch(snes_, &linesearch); petsc_err(perr);
        #endif
        perr = SNESLineSearchSetType(linesearch, "basic"); petsc_err(perr);
        perr = SNESLineSearchSetComputeNorms(linesearch, PETSC_FALSE);
        #else
        perr = SNESLineSearchSet(snes_, SNESLineSearchNoNorms, PETSC_NULL); petsc_err(perr); 
        #endif
      }
      else
      {
        tf_err("Unknown snes ls type.", "Requested ls type: %s", lstype.c_str());
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
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR != 2
        tf_err("Cannot set snes ls min_lambda with PETSc != 3.2.", "PETSc version: %d.%d",
               PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
        #endif
      }
      double minlambda;
      serr = Spud::get_option(buffer.str(), minlambda, 1.e-12);
      spud_err(buffer.str(), serr);

      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #if PETSC_VERSION_MINOR > 2
      SNESLineSearch linesearch;
      #if PETSC_VERSION_MINOR > 3
      perr = SNESGetLineSearch(snes_, &linesearch); petsc_err(perr);
      #else
      perr = SNESGetSNESLineSearch(snes_, &linesearch); petsc_err(perr);
      #endif
      if (lstype == "cubic" || lstype == "quadratic")
      {
        perr = SNESLineSearchBTSetAlpha(linesearch, alpha); petsc_err(perr);// FIXME: assumes using bt
      }
      perr = SNESLineSearchSetTolerances(linesearch, PETSC_DEFAULT, maxstep, PETSC_DEFAULT,
                                        PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
      petsc_err(perr);
      #else
      perr = SNESLineSearchSetParams(snes_, alpha, maxstep, minlambda); petsc_err(perr);
      #endif
      #else
      perr = SNESLineSearchSetParams(snes_, alpha, maxstep); petsc_err(perr);
      #endif
       
    }
    else
    {
      perr = SNESSetType(snes_, snestype.c_str()); petsc_err(perr); 
      perr = SNESSetFromOptions(snes_); petsc_err(perr);               // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)
    }

    buffer.str(""); buffer << optionpath() << "/type/snes_type/convergence_test/name";
    if (Spud::have_option(buffer.str()))
    {
      std::string convtest;
      serr = Spud::get_option(buffer.str(), convtest);
      if (convtest=="skip")
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
        perr = SNESSetConvergenceTest(snes_, SNESSkipConverged, 
                                            PETSC_NULL, PETSC_NULL);
        #else
        perr = SNESSetConvergenceTest(snes_, SNESConvergedSkip, 
                                            PETSC_NULL, PETSC_NULL);
        #endif
      }
    }

    buffer.str(""); buffer << optionpath() 
                                        << "/type/monitors/residual";
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 7
      perr = SNESMonitorSet(snes_, SNESMonitorDefault,               // set a snes residual monitor
                                            PETSC_NULL, PETSC_NULL); 
      #else
      PetscViewerAndFormat *vf;
      PetscViewer       viewer;
      perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)snes_),&viewer);
      petsc_err(perr);
      perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf);
      petsc_err(perr);
      perr = SNESMonitorSet(snes_, (PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))SNESMonitorDefault,               // set a snes residual monitor
                                            vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); 
      #endif
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath() 
                                  << "/type/monitors/solution_graph";
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 7
      perr = SNESMonitorSet(snes_, SNESMonitorSolution,              // set a snes solution monitor (graph)
                                            PETSC_NULL, PETSC_NULL); 
      #else
      PetscViewerAndFormat *vf;
      PetscViewer       viewer;
      perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)snes_),&viewer);
      petsc_err(perr);
      perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf);
      petsc_err(perr);
      perr = SNESMonitorSet(snes_, (PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))SNESMonitorSolution,               // set a snes residual monitor
                                            vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); 
      #endif
      petsc_err(perr);
    }

    if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
    {
      snesmctx_.solver = this;
      if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
      {
        buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                               << (*system()).name() << "_" 
                               << name() << "_snes.conv";
        convfile_.reset( new ConvergenceFile(buffer.str(),
                                      (*(*system_).mesh()).mpi_comm(),// allocate the file but don't write the header yet as the
                                      &(*(*system()).bucket()), (*system()).name(), name()) ); // bucket isn't complete
      }
      perr = SNESMonitorSet(snes_, SNESCustomMonitor,                // set a custom snes monitor
                                            &snesmctx_, PETSC_NULL); 
      petsc_err(perr);
    }

    perr = SNESSetTolerances(snes_, atol_, rtol_, stol_, maxits_,    // from the data we collected in the base data fill set-up the
                                                            maxfes_);// snes tolerances
    petsc_err(perr);

    perr = SNESGetKSP(snes_, &ksp_); petsc_err(perr);                  // we always have at least one ksp so use the solverbucket ksp
                                                                     // to start setting up the ksp inside the snes
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // the ksp solver path
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // can then be used to fill the ksp data

    buffer.str(""); buffer << optionpath() << "/type/monitors/view_snes";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESView(snes_, 
             PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm())); 
      petsc_err(perr);                                                 // turn on snesview so we get some debugging info
    }

  }
  else if (type()=="Picard")                                         // if this is a picard solver
  {

    perr = KSPCreate((*(*system_).mesh()).mpi_comm(), &ksp_); 
    petsc_err(perr);                                                   // create a ksp object from the variable in the solverbucket

    if (Spud::have_option(optionpath()+"/type/monitors/convergence_file"))
    {
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_" 
                             << name() << "_picard.conv";
      convfile_.reset( new ConvergenceFile(buffer.str(), 
                                    (*(*system_).mesh()).mpi_comm(), // allocate the file but don't write the header yet as the
                                    &(*(*system()).bucket()), (*system()).name(), name()) );   // bucket isn't complete
    }

    if (bilinearpc_)
    {                                                                // if there's a pc associated
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = KSPSetOperators(ksp_, (*matrix_).mat(),                // set the ksp operators with two matrices
                                   (*matrixpc_).mat(), 
                                   SAME_NONZERO_PATTERN); 
      #else
      perr = KSPSetOperators(ksp_, (*matrix_).mat(),                // set the ksp operators with two matrices
                                   (*matrixpc_).mat()); 
      #endif
      petsc_err(perr);
    }
    else
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = KSPSetOperators(ksp_, (*matrix_).mat(),                // set the ksp operators with the same matrices
                                   (*matrix_).mat(), 
                                   SAME_NONZERO_PATTERN); 
      #else
      perr = KSPSetOperators(ksp_, (*matrix_).mat(),                // set the ksp operators with the same matrices
                                   (*matrix_).mat()); 
      #endif
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // figure out the linear solver optionspath
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // fill the ksp data

    buffer.str(""); buffer << optionpath() << "/type/linear_solver/monitors/view_ksp";
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPView(ksp_, 
             PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm())); 
      petsc_err(perr);                                                 // turn on kspview so we get some debugging info
    }

  }
  else                                                               // don't know how we got here
  {
    tf_err("Unknown solver type.", "Solver type: %s", type().c_str());
  }

  create_nullspace();                                                // this should be safe to call now as null spaces will
                                                                     // have been initialized from all the solvers

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps (with an optionpath as well)
//*******************************************************************|************************************************************//
void SpudSolverBucket::register_form(Form_ptr form, 
                                      const std::string &name, 
                                      std::string optionpath)
{
  Form_hash_it f_it = forms_.get<om_key_hash>().find(name);                                  // check if name exists
  if (f_it != forms_.get<om_key_hash>().end())
  {
    tf_err("Form already exists in solver.", "Form name: %s. Solver name: %s.", name.c_str(), name_.c_str());
  }
  else
  {
    forms_.insert(om_item<const std::string, Form_ptr>(name, form));                                  // if not, insert form pointer into data map
    form_optionpaths_.insert(om_item<const std::string, std::string>(name, optionpath));                            // and do the same for its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string containing the named form's optionpath
//*******************************************************************|************************************************************//
const std::string SpudSolverBucket::fetch_form_optionpath(const std::string &name) const
{
  string_hash_it s_it = form_optionpaths_.get<om_key_hash>().find(name);
                                                                     // check if the name already exists
  if (s_it == form_optionpaths_.get<om_key_hash>().end())
  {
    tf_err("Form does not exist in solver.", "Form name: %s. Solver name: %s.", name.c_str(), name_.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the form_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudSolverBucket::form_optionpaths_begin()
{
  return form_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the form_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudSolverBucket::form_optionpaths_begin() const
{
  return form_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the form_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudSolverBucket::form_optionpaths_end()
{
  return form_optionpaths_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the form_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudSolverBucket::form_optionpaths_end() const
{
  return form_optionpaths_.get<om_key_seq>().end();
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

  for ( string_const_it s_it = form_optionpaths_begin(); 
                            s_it != form_optionpaths_end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << 
                                (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// checkpoint the options file
//*******************************************************************|************************************************************//
void SpudSolverBucket::checkpoint_options_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  if (solve_location() == SOLVE_START)                               // do not solve again if this system was only meant to be
  {                                                                  // solved once
    std::string location = "never";
    buffer.str(""); buffer << optionpath() << "/solve/name";
    serr = Spud::set_option_attribute(buffer.str(), location);
    spud_err(buffer.str(), serr);
  }

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

  buffer.str(""); buffer << optionpath() << 
                                "/type/monitors/norms";
  monitornorms_ = Spud::have_option(buffer.str());

  std::string location;
  buffer.str(""); buffer << optionpath() << "/solve/name";
  serr = Spud::get_option(buffer.str(), location);
  spud_err(buffer.str(), serr);
  if (location=="in_timeloop")
  {
    solve_location_ = SOLVE_TIMELOOP;
  }
  else if (location=="at_start")
  {
    solve_location_ = SOLVE_START;
  }
  else if (location=="with_diagnostics")
  {
    solve_location_ = SOLVE_DIAGNOSTICS;
  }
  else if (location=="never")
  {
    solve_location_ = SOLVE_NEVER;
  }
  else
  {
    tf_err("Unknown solve location.", "Solver: %s::%s", (*system()).name().c_str(), name_.c_str());
  }

  solved_.reset( new bool(false) );                                  // assume the solver hasn't been solved yet

  sp_   = PETSC_NULL;                                                // initialize in case we don't get a chance
  ksp_  = PETSC_NULL;                                                // to do this later
  snes_ = PETSC_NULL;

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
    residual_ = linear_;                                             // the linear form is the residual for SNES

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
    tf_err("Unknown solver type.", "Solver type: %s", type().c_str());
  }

  for (Form_const_it f_it = forms_begin(); f_it != forms_end(); f_it++)
  {
    if(boost::algorithm::ends_with((*f_it).first, "SchurPC"))
    {
      register_solverform((*f_it).second, (*f_it).first);
      solverident_zeros_[(*f_it).first] = Spud::have_option(fetch_form_optionpath((*f_it).first)+"/ident_zeros");
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
  work_.reset( new dolfin::PETScVector(*std::dynamic_pointer_cast<dolfin::PETScVector>((*(*system_).function()).vector())) ); 
  (*work_).zero();

  dolfin::SystemAssembler sysassembler(bilinear_, linear_, 
                                    (*system_).bcs());
  sysassembler.keep_diagonal = true;
  matrix_.reset(new dolfin::PETScMatrix);                            // allocate the matrix
  sysassembler.assemble(*matrix_);

  if(bilinearpc_)                                                    // do we have a pc form?
  {
    dolfin::SystemAssembler sysassemblerpc(bilinearpc_, linear_,
                                           (*system_).bcs());
    sysassemblerpc.keep_diagonal = true;
    matrixpc_.reset(new dolfin::PETScMatrix);                        // allocate the matrix
    sysassemblerpc.assemble(*matrixpc_);
  }

  for (Form_const_it f_it = solverforms_begin(); 
                     f_it != solverforms_end(); f_it++)
  {
    dolfin::SystemAssembler sysassemblerform((*f_it).second, linear_,
                                             (*system_).bcs());
    sysassemblerform.keep_diagonal = true;
    PETScMatrix_ptr solvermatrix;
    solvermatrix.reset(new dolfin::PETScMatrix);
    sysassemblerform.assemble(*solvermatrix);
    solvermatrices_[(*f_it).first] = solvermatrix;
  }

  dolfin::Assembler assembler;
   
  rhs_.reset(new dolfin::PETScVector);                               // allocate the rhs
  assembler.assemble(*rhs_, *linear_);

  res_.reset(new dolfin::PETScVector);                               // allocate the residual
  assembler.assemble(*res_, *residual_);

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
  petsc_err(perr);

  std::string iterative_method;
  buffer.str(""); buffer << optionpath << "/iterative_method/name";  // iterative method (gmres, fgmres, cg etc.)
  serr = Spud::get_option(buffer.str(), iterative_method); 
  spud_err(buffer.str(), serr);

  perr = KSPSetType(ksp, iterative_method.c_str()); petsc_err(perr);   // set the ksp type to the iterative method

  perr = KSPSetFromOptions(ksp);                                     // do this now so that options will be overwritten by options
                                                                     // file

  PC pc;
  perr = KSPGetPC(ksp, &pc); petsc_err(perr);                          // get the pc from the ksp

  uint parent_offset = 0;
  if (parent_indices)
  {
    parent_offset = dolfin::MPI::global_offset((*(*system_).mesh()).mpi_comm(),
                                               (*parent_indices).size(), true);
  }

  fill_pc_(optionpath, pc,
           prefix,
           parent_offset, parent_indices);

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
      petsc_err(perr);
    }
    else
    {
      perr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); 
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/restart";
    if (Spud::have_option(buffer.str()))
    {
      PetscInt restart;
      serr = Spud::get_option(buffer.str(), restart);
      spud_err(buffer.str(), serr);
      perr = KSPGMRESSetRestart(ksp, restart);
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath << 
                "/iterative_method/monitors/preconditioned_residual";// monitor the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 7
      perr = KSPMonitorSet(ksp, KSPMonitorDefault, 
                                            PETSC_NULL, PETSC_NULL); 
      #else
      PetscViewerAndFormat *vf;
      PetscViewer       viewer;
      perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ksp),&viewer);
      petsc_err(perr);
      perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf);
      petsc_err(perr);
      perr = KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorDefault,
                                            vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); 
      #endif
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath << 
                          "/iterative_method/monitors/true_residual";// monitor the true residual (more expensive)
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 7
      perr = KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, 
                                            PETSC_NULL, PETSC_NULL); 
      #else
      PetscViewerAndFormat *vf;
      PetscViewer       viewer;
      perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ksp),&viewer);
      petsc_err(perr);
      perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf);
      petsc_err(perr);
      perr = KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorTrueResidualNorm,
                                            vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); 
      #endif
      petsc_err(perr);
    }

    buffer.str(""); buffer << optionpath << 
          "/iterative_method/monitors/preconditioned_residual_graph";// plot a graph of the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 4
      tf_err("Preconditioned residual graph not available", "Not supported with PETSc > 3.4.");
      #elif PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
      perr = KSPMonitorSet(ksp, KSPMonitorLGResidualNorm, 
                                             PETSC_NULL, PETSC_NULL); 
      #else
      perr = KSPMonitorSet(ksp, KSPMonitorLG, 
                                             PETSC_NULL, PETSC_NULL); 
      #endif
      petsc_err(perr);
    }

    if (Spud::have_option(optionpath+"/iterative_method/monitors/convergence_file"))
    {
      kspmctx_.solver = this;
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_" 
                             << name() << "_ksp.conv";
      kspconvfile_.reset( new KSPConvergenceFile(buffer.str(), 
                                    (*(*system_).mesh()).mpi_comm(), // allocate the file but don't write the header yet as the
                                    &(*(*system()).bucket()), (*system()).name(), name()) );   // bucket isn't complete
      perr = KSPMonitorSet(ksp, KSPCustomMonitor, 
                                             &kspmctx_, PETSC_NULL); 
      petsc_err(perr);
    }

    if (Spud::have_option(optionpath+"/iterative_method/monitors/test_null_space"))
    {
      kspmctx_.solver = this;
      perr = KSPMonitorSet(ksp, KSPNullSpaceMonitor, 
                                             &kspmctx_, PETSC_NULL); 
      petsc_err(perr);
    }

    perr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  }

  buffer.str(""); buffer << optionpath << "/remove_null_space";      // removing a (or multiple) null space(s)
  if (Spud::have_option(buffer.str()))
  {
    MatNullSpace SP;                                                 // create a set of nullspaces in a null space object
    fill_nullspace_(buffer.str(), SP, parent_offset, parent_indices);

    Mat mat;
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = KSPGetOperators(ksp, &mat, PETSC_NULL, PETSC_NULL);
    #else
    perr = KSPGetOperators(ksp, &mat, PETSC_NULL);
    #endif
    petsc_err(perr);

    perr = MatSetNullSpace(mat, SP); petsc_err(perr);

    #if PETSC_VERSION_MAJOR < 4 && PETSC_VERSION_MINOR < 6
    perr = KSPSetNullSpace(ksp, SP); petsc_err(perr);                  // attach it to the ksp
    #endif

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatNullSpaceDestroy(&SP); petsc_err(perr);                  // destroy the null space object, necessary?
    #else
    perr = MatNullSpaceDestroy(SP); petsc_err(perr);                   // destroy the null space object, necessary?
    #endif
  }

}

//*******************************************************************|************************************************************//
// fill a pc object from the options tree (this routine may be called recursively for ksp, bjacobi, asm and fieldsplit pc types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_pc_(const std::string &optionpath, PC &pc, 
                                const std::string prefix,
                                const uint &parent_offset,
                                const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  perr = PCSetOptionsPrefix(pc, prefix.c_str());                     // set the pc options prefix
  petsc_err(perr);

  std::string preconditioner;
  buffer.str(""); buffer << optionpath << "/preconditioner/name";    // preconditioner type
  serr = Spud::get_option(buffer.str(), preconditioner);
  spud_err(buffer.str(), serr);

  perr = PCSetType(pc, preconditioner.c_str()); petsc_err(perr);       // set its type (read from options earlier)

  perr = PCSetFromOptions(pc); petsc_err(perr);                        // do this now so they can be overwritten

  if (preconditioner=="ksp")                                         // if the pc is itself a ksp
  {
    buffer.str(""); buffer << optionpath << 
                                    "/preconditioner/linear_solver";
    KSP subksp;                                                      // create a subksp from this pc
    perr = PCKSPGetKSP(pc, &subksp); petsc_err(perr);
    fill_ksp_(buffer.str(), subksp, prefix, parent_indices);         // recursively fill the ksp data (i.e. go back to this routine)
  }
  else if (preconditioner=="lsc")                                    // would be nice to put subksp options for lsc here
  {
    perr = PCSetUp(pc); petsc_err(perr);                               // this is just necessary in case we view it later
  }
  else if (preconditioner=="fieldsplit")                             // if the pc is a fieldsplit
  {
    buffer.str(""); buffer << optionpath << "/preconditioner";
    fill_pc_fieldsplit_(buffer.str(), pc, prefix, parent_offset, 
                                                  parent_indices);   // fill the fieldsplit data (will end up back here again
                                                                     // eventually)
  }
  else if (preconditioner=="lu")                                     // if the pc is direct
  {
    std::string factorization_package;                               // we get to choose a factorization package
    buffer.str(""); buffer << optionpath << 
                        "/preconditioner/factorization_package/name";
    serr = Spud::get_option(buffer.str(), factorization_package); 
    spud_err(buffer.str(), serr);

    perr = PCFactorSetMatSolverType(pc, 
                                      factorization_package.c_str()); 
    petsc_err(perr);

  }
  else if (preconditioner=="hypre")
  {
    std::string hypre_type;                                          // we get to choose the hypre type
    buffer.str(""); buffer << optionpath << 
                        "/preconditioner/hypre_type/name";
    serr = Spud::get_option(buffer.str(), hypre_type); 
    spud_err(buffer.str(), serr);

#ifdef PETSC_HAVE_HYPRE
    perr = PCHYPRESetType(pc, hypre_type.c_str()); petsc_err(perr);
#else
    tf_err("Hypre not available", "PETSc not compiled with hypre.");
#endif
  }
  else if ((preconditioner=="bjacobi")||(preconditioner=="asm"))
  {
    perr = PCSetUp(pc); petsc_err(perr);                               // call this before subpc can be retrieved

    KSP *subksp;
    if(preconditioner=="bjacobi")
    {
      perr = PCBJacobiGetSubKSP(pc, PETSC_NULL, PETSC_NULL, &subksp);
      petsc_err(perr);
    }
    else if (preconditioner=="asm")
    {
      perr = PCASMGetSubKSP(pc, PETSC_NULL, PETSC_NULL, &subksp);
      petsc_err(perr);
    }
    else
    {
      tf_err("Unknown preconditioner.", "Preconditioner: %s", preconditioner.c_str());
    }

    PC subpc;
    perr = KSPGetPC(*subksp, &subpc); petsc_err(perr);                  // get the sub pc from the sub ksp
    fill_pc_(optionpath+"/preconditioner", subpc,
             prefix, 
             parent_offset, parent_indices);

  }

  buffer.str(""); buffer << optionpath 
                                << "/preconditioner/near_null_space";// set a (or multiple) near null space(s)
  if (Spud::have_option(buffer.str()))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
    Mat pmat;
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = PCGetOperators(pc, PETSC_NULL, &pmat, PETSC_NULL);
    #else
    perr = PCGetOperators(pc, PETSC_NULL, &pmat);
    #endif
    petsc_err(perr);

    MatNullSpace SP;                                                 // create a set of nullspaces in a null space object
    fill_nullspace_(buffer.str(), SP, parent_offset, parent_indices);

    perr = MatSetNearNullSpace(pmat, SP); petsc_err(perr);

    perr = MatNullSpaceDestroy(&SP); petsc_err(perr);                  // destroy the null space object, necessary?
    #else
    tf_err("Can only set near null spaces with petsc > 3.2.", "PETSc version: %d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
    #endif
  }

  if ((preconditioner=="ml")||(preconditioner=="gamg"))
  {
    perr = PCSetUp(pc); petsc_err(perr);                               // need to call this to prevent seg fault on kspview
  }                                                                  // BUT it has to happen after the near null space is set

}

//*******************************************************************|************************************************************//
// fill a pc object from the options tree assuming its a fieldsplit pc
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_pc_fieldsplit_(const std::string &optionpath, 
                                           PC &pc, const std::string prefix,
                                           const uint &parent_offset,
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
    if (i==0)
    {
      fill_indices_values_by_field_(buffer.str(),                                // setup an IS for each fieldsplit
                                    indices, NULL, 
                                    parent_indices, 
                                    NULL);

      prev_indices = indices;                                        // should already be sorted so don't bother doing it again
    }
    else
    {
      fill_indices_values_by_field_(buffer.str(),
                                    indices, NULL, 
                                    parent_indices,
                                    &prev_indices);

      prev_indices.insert(prev_indices.end(), indices.begin(), indices.end());
      std::sort(prev_indices.begin(), prev_indices.end());           // sort the vector of prev_indices
    }

    IS is = convert_vector_to_is((*(*system_).mesh()).mpi_comm(), 
                                 indices,
                                 parent_offset, parent_indices);

    buffer.str(""); buffer << optionpath << 
                                        "/fieldsplit[" 
                                         << i << "]/name";
    std::string fsname;
    serr = Spud::get_option(buffer.str(), fsname);
    spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << "/fieldsplit[" 
                                    << i << "]/monitors/view_index_set";
    if (Spud::have_option(buffer.str()))
    {
      if (dolfin::MPI::rank((*(*system_).mesh()).mpi_comm())==0)
      {
        log(INFO, "ISView: %s (%s)", 
                                    fsname.c_str(), optionpath.c_str());
      }
      perr = ISView(is, 
                PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm()));
      petsc_err(perr);                                                   // isview?
    }

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    buffer.str(""); buffer << prefix << fsname;
    perr = PCFieldSplitSetIS(pc, buffer.str().c_str(), is);          // set the fs using that IS
    petsc_err(perr);
    perr = ISDestroy(&is); petsc_err(perr);                            // destroy the IS, necessary?
    #else
    perr = PCFieldSplitSetIS(pc, is); petsc_err(perr);                 // set the fs using that IS
    perr = ISDestroy(is); petsc_err(perr);                             // destroy the IS, necessary?
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
    petsc_err(perr);
  }
  else if (ctype == "multiplicative")                                // multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE); 
    petsc_err(perr);
  }
  else if (ctype == "symmetric_multiplicative")                      // symmetric multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, 
                          PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE); 
    petsc_err(perr);
  }
  else if (ctype == "special")                                       // special fieldsplit (whatever that means!)
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SPECIAL); 
    petsc_err(perr);
  }
  else if (ctype == "schur")                                         // schur fieldsplit
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR); 
    petsc_err(perr);
    
    std::string ftype;
    buffer.str(""); buffer << optionpath << 
                    "/composite_type::schur/factorization_type/name";// schur factorization type
    serr = Spud::get_option(buffer.str(), ftype);                    // sadly no string based interface provided to this (I think)
    spud_err(buffer.str(), serr);                                    // so hard code an if block
    if (ftype == "full")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
      perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL); 
      petsc_err(perr);
      #endif
    }
    else
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 2
      if (ftype == "upper")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER); 
        petsc_err(perr);
      }
      else if (ftype == "lower")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER); 
        petsc_err(perr);
      }
      else if (ftype == "diag")
      {
        perr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG); 
        petsc_err(perr);
      }
      else
      {
        tf_err("Unknown PCFieldSplitSchurFactType.", "PCFieldSplitSchutFactType: %s", ftype.c_str());
      }
      #else
      tf_err("Can only set schur factorization_type to anything other than full with petsc > 3.2.", 
             "PETSc version: %d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
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
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_A11, PETSC_NULL);
      #else
      perr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_A11, PETSC_NULL);
      #endif
      #else
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_DIAG, PETSC_NULL);
      #else
      perr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_DIAG, PETSC_NULL);
      #endif
      #endif
      petsc_err(perr);
    }
    else if (ptype == "self")
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_SELF, PETSC_NULL);
      #else
      perr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELF, PETSC_NULL);
      #endif
      petsc_err(perr);
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

      IS is;
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      perr = ISCreateGeneral((*(*system_).mesh()).mpi_comm(), n, 
                            petscindices, PETSC_OWN_POINTER, &is);// create the general index set based on the indices
      #else
      perr = ISCreateGeneral((*(*system_).mesh()).mpi_comm(), n, 
                                               petscindices, &is);// create the general index set based on the indices
      #endif
      petsc_err(perr);
      
      buffer.str(""); buffer << optionpath << 
                 "/composite_type::schur/schur_preconditioner::user";
      if (Spud::have_option(buffer.str()+"/monitors/view_index_set"))
      {
        if (dolfin::MPI::rank((*(*system_).mesh()).mpi_comm())==0)
        {
          std::string isname = prefix+"SchurPC";
          log(INFO, "ISView: %s (%s)", 
                               isname.c_str(), buffer.str().c_str());
        }
        perr = ISView(is, 
              PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm())); 
        petsc_err(perr);                                               // isview?
      }

      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      #else
      PetscFree(petscindices);                                       // free the PetscInt array of indices
      #endif

      solverindexsets_[prefix+"SchurPC"] = is;

      Mat submatrix;
      perr = MatCreateSubMatrix((*solvermatrices_[prefix+"SchurPC"]).mat(), is, is, MAT_INITIAL_MATRIX, &submatrix);
      petsc_err(perr);

      solversubmatrices_[prefix+"SchurPC"] = submatrix;
      
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      perr = PCFieldSplitSchurPrecondition(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, submatrix);
      #else
      perr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, submatrix);
      #endif
      petsc_err(perr);
    }
    else
    {
      tf_err("Unknown PCFieldSplitSchurPreconditionType.", 
             "PCFieldSplitSchurPreconditionType: %s", ptype.c_str());
    }
    #else
    tf_err("Schur fieldsplits only supprted with petsc > 3.1.", 
           "PETSc version: %d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
    #endif
  }
  else                                                               // unknown (to buckettools) fieldsplit composite type
  {
    tf_err("Unknown PCCompositeType.", 
           "PCCompositeType: %s", ctype.c_str());
  }
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
  perr = PCSetUp(pc); petsc_err(perr);                                 // call this before subksp can be retrieved
                                                                     // this is only necessary for schur fieldsplits and
                                                                     // as we're not supporting these with petsc 3.1 we don't
                                                                     // bother otherwise
  #endif

  KSP *subksps;                                                      // setup the fieldsplit subksps
  PetscInt nsubksps;
  perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); 
  petsc_err(perr); 

  assert(nsubksps==nsplits);

  for (uint i = 0; i < nsplits; i++)                                 // loop over the splits again
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = KSPSetNormType(subksps[i], KSP_NORM_DEFAULT);             // these two calls are necessary because the normtype
    petsc_err(perr);                                                   // and pc_side will have been set inappropriately in the 
    perr = KSPSetPCSide(subksps[i], PC_SIDE_DEFAULT);                // above call to PCSetUp.  now we undo this damage by resetting
    petsc_err(perr);                                                   // them to default (since the schema doesn't allow them to be
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
// Fill a petsc IS object from fields
// IS's may be set up by field name, components of the field, regions of the domain of the field  and surfaces of the domain of the
// field.
// Additionally the resulting IS is checked for consistency with any parents or siblings in the tree.
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_indices_values_by_field_(const std::string &optionpath,
                                                     std::vector<uint> &child_indices, 
                                                     PETScVector_ptr values,
                                                     const std::vector<uint>* parent_indices,
                                                     const std::vector<uint>* sibling_indices)
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  buffer.str(""); buffer << optionpath << "/field";                  // loop over the fields used to describe this IS
  int nfields = Spud::option_count(buffer.str());

  PETScVector_ptr tmp_values;
  if (values)
  {
    tmp_values.reset( new dolfin::PETScVector(*values) );
  }
 
  child_indices.clear();
  if (nfields==0)                                                    // if no fields have been specified... **no fields**
  {
    std::pair<uint, uint> ownership_range =                          // the parallel ownership range of the system functionspace
            (*(*(*system_).functionspace()).dofmap()).ownership_range();
    for (uint i = ownership_range.first; i < ownership_range.second; i++)
    {
      child_indices.push_back(i);
      if (tmp_values)
      {
        (*tmp_values).setitem(i, 1.0);
      }
    }
  }
  else
  {

    const bool mixedsystem = ((*system_).fields_size()>1);

    for (uint i = 0; i < nfields; i++)                               // loop over the fields that have been specified
    {

      std::string fieldname;
      buffer.str(""); buffer << optionpath << "/field[" << i <<      // get the field name
                                                            "]/name";
      serr = Spud::get_option(buffer.str(), fieldname); 
      spud_err(buffer.str(), serr);
      
      FunctionBucket_ptr field = (*system_).fetch_field(fieldname);  // using the name, get the
      const int fieldindex = (*field).index();                       // field index
      const std::vector<std::size_t> fieldshape = (*field).shape();  // field shape
      const bool fieldsymmetric = (*field).symmetric();              // field symmetric

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
      dolfin::Expression* value_exp = NULL;
      double *value_const = NULL;
      buffer.str(""); buffer << optionpath << "/field[" << i << "]";
      field_restrictions_(buffer.str(),                              // get any restrictions on the IS
                          components, region_ids, boundary_ids, 
                          value_exp, value_const, 
                          fieldshape, fieldsymmetric);
     
      std::vector<uint> tmp_indices;
      tmp_indices  = functionspace_dofs_values(functionspace, 
                                               (*system_).celldomains(),
                                               (*system_).facetdomains(),
                                               components, region_ids, 
                                               boundary_ids,
                                               tmp_values, value_exp, value_const);

      child_indices.insert(child_indices.end(), tmp_indices.begin(), tmp_indices.end());

      destroy_field_restrictions_(components, region_ids, 
                                  boundary_ids, value_exp, 
                                  value_const);
    }
  }

  restrict_indices(child_indices, 
                   (*system_).functionspace(), 
                   parent_indices, 
                   sibling_indices);

  if(values)
  {
    restrict_values(values, tmp_values, child_indices);
  }

}

//*******************************************************************|************************************************************//
// fill a petsc nullspace object using the options in the optionpath provided
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_nullspace_(const std::string &optionpath, MatNullSpace &SP,
                                       const uint &parent_offset, const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  buffer.str(""); buffer << optionpath << "/null_space";
  int nnulls = Spud::option_count(buffer.str());                     // how many null spaces?
  buffer.str(""); buffer << optionpath <<  "/null_space/remove_from_rhs/only_remove_from_rhs";
  int rhsnnulls = Spud::option_count(buffer.str());

  PETScVector_ptr sysvec = std::dynamic_pointer_cast<dolfin::PETScVector>((*(*system_).function()).vector());

  Vec vecs[nnulls-rhsnnulls];                                        // and here (for the petsc interface)

  for (uint i = 0; i<nnulls; i++)                                    // loop over the nullspaces
  {

    std::vector<uint> indices;                                       // the indices of the system vector used in the null space
    PETScVector_ptr nullvec( new dolfin::PETScVector(*sysvec) );
    (*nullvec).zero();                                               // create a null vector (with the size and attributes of the
                                                                     // system)

    buffer.str(""); buffer << optionpath <<                          // optionpath of the nullspace
                      "/null_space[" << i << "]";
    fill_indices_values_by_field_(buffer.str(), indices, nullvec,    // create a vector describing the nullspace based on this optionpath 
                                  parent_indices);                   // (no siblings as null spaces can overlap)

    perr = VecNormalize((*nullvec).vec(), PETSC_NULL); petsc_err(perr);// normalize the null space vector

    buffer.str(""); buffer << optionpath << "/null_space[" 
                                 << i << "]/remove_from_rhs";
    if (Spud::have_option(buffer.str()))
    {
      nullspacevectors_.push_back(nullvec);                          // permanently store system wide null vector
    }

    IS sis = convert_vector_to_is((*(*system_).mesh()).mpi_comm(),   // this IS indexes from the system vector into the child
                                 indices);

    IS pis = convert_vector_to_is((*(*system_).mesh()).mpi_comm(),   // this IS indexes from the parent vector into the child
                                 indices,
                                 parent_offset, parent_indices);

    buffer.str(""); buffer << optionpath << "/null_space[" 
                                 << i << "]/monitors/view_index_set";
    if (Spud::have_option(buffer.str()))
    {
      buffer.str(""); buffer << optionpath << 
                                          "/null_space[" 
                                           << i << "]/name";
      std::string nsname;
      serr = Spud::get_option(buffer.str(), nsname);
      spud_err(buffer.str(), serr);

      if (dolfin::MPI::rank((*(*system_).mesh()).mpi_comm())==0)
      {
        log(INFO, "ISView: nullspace (system) %s (%s)", 
                                    nsname.c_str(), optionpath.c_str());
      }
      perr = ISView(sis, 
                PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm()));
      petsc_err(perr);                                                 // isview?

      if (dolfin::MPI::rank((*(*system_).mesh()).mpi_comm())==0)
      {
        log(INFO, "ISView: nullspace (subsystem) %s (%s)", 
                                    nsname.c_str(), optionpath.c_str());
      }
      perr = ISView(pis, 
                PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm()));
      petsc_err(perr);                                                 // isview?
    }

    buffer.str(""); buffer << optionpath << "/null_space[" 
                                 << i << "]/remove_from_rhs/only_remove_from_rhs";
    if (Spud::have_option(buffer.str()))
    {
      continue;
    }

    perr = VecCreate((*(*system_).mesh()).mpi_comm(), &vecs[i]);
    petsc_err(perr);
    if (parent_indices)
    {
      perr = VecSetSizes(vecs[i], (*parent_indices).size(), 
                                             PETSC_DECIDE);
      petsc_err(perr);
    }
    else
    {
      perr = VecSetSizes(vecs[i], (*nullvec).local_size(), 
                                             PETSC_DECIDE);
      petsc_err(perr);
    }
    perr = VecSetUp(vecs[i]);
    perr = VecSet(vecs[i], (double) 0.0); petsc_err(perr);

    VecScatter scatter;
    perr = VecScatterCreate((*nullvec).vec(), sis, 
                            vecs[i], pis, 
                            &scatter);
    petsc_err(perr);
    perr = VecScatterBegin(scatter, 
                           (*nullvec).vec(), vecs[i], 
                           INSERT_VALUES, SCATTER_FORWARD);
    petsc_err(perr);
    perr = VecScatterEnd(scatter,
                         (*nullvec).vec(), vecs[i],
                         INSERT_VALUES, SCATTER_FORWARD);
    petsc_err(perr);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1          // necessary or taken care of when object leaves scope?
    perr = VecScatterDestroy(&scatter); petsc_err(perr);      
    perr = ISDestroy(&sis); petsc_err(perr);
    perr = ISDestroy(&pis); petsc_err(perr);
    #else
    perr = VecScatterDestroy(scatter); petsc_err(perr);
    perr = ISDestroy(sis); petsc_err(perr);
    perr = ISDestroy(pis); petsc_err(perr);
    #endif

  }

  orthonormalize_petsc_vecs_(vecs, nnulls-rhsnnulls);

  perr = MatNullSpaceCreate((*(*system_).mesh()).mpi_comm(), 
                            PETSC_FALSE, nnulls-rhsnnulls, vecs, &SP); 
  petsc_err(perr);

  buffer.str(""); buffer << optionpath << 
                    "/monitors/view_null_space";                     // view the null space for debugging
  if (Spud::have_option(buffer.str()))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatNullSpaceView(SP, 
           PETSC_VIEWER_STDOUT_((*(*system_).mesh()).mpi_comm())); 
    petsc_err(perr);
    #else
    log(WARNING, "Cannot set view_null_space monitor with PETSc < 3.2.");
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
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
  fill_bound_(buffer.str(), ub, SNES_VI_INF);
  #else
  fill_bound_(buffer.str(), ub, PETSC_INFINITY);
  #endif

  PETScVector_ptr lb;
  buffer.str(""); buffer << optionpath() << "/type/snes_type/constraints/lower_bound";
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
  fill_bound_(buffer.str(), lb, SNES_VI_NINF);
  #else
  fill_bound_(buffer.str(), lb, PETSC_NINFINITY);
  #endif

  perr = SNESVISetVariableBounds(snes_, (*lb).vec(), (*ub).vec());
  petsc_err(perr);
                                                                     // FIXME: UGLY HACK: our constant bounds will be overwritten by the dm
                                                                     // in SNESSetUp so let's stop it from doing that by attaching a
                                                                     // dummy variable bounds computation - does nothing!
  perr = SNESVISetComputeVariableBounds(snes_, SNESVIDummyComputeVariableBounds);
  petsc_err(perr);
  #else
  tf_err("SNES VI only available with petsc > 3.1.", 
         "PETSc version: %d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);
  #endif

}

//*******************************************************************|************************************************************//
// fill petsc vectors describing the bounds on fields using the optionpath provided
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_bound_(const std::string &optionpath, PETScVector_ptr &bound, const double &background_value)
{

  PETScVector_ptr sysvec = std::dynamic_pointer_cast<dolfin::PETScVector>((*(*system_).function()).vector());
  bound.reset( new dolfin::PETScVector(*sysvec) );

  std::vector<double> background((*bound).local_size(), background_value);
  (*bound).set_local(background);
  (*bound).apply("insert");
  background.clear();

  if (Spud::have_option(optionpath))
  {
    std::vector<uint> indices;
    fill_indices_values_by_field_(optionpath, indices, bound);
    indices.clear();
  }

}

//*******************************************************************|************************************************************//
// using the optionpath set up any restrictions we place on the field IS (i.e. region ids, components or boundary_ids)
//*******************************************************************|************************************************************//
void SpudSolverBucket::field_restrictions_(const std::string &optionpath,
                                           std::vector<int>* &components,
                                           std::vector<int>* &region_ids,
                                           std::vector<int>* &boundary_ids,
                                           dolfin::Expression* &value_exp,
                                           double* &value_const,
                                           const std::vector<std::size_t> &fieldshape,
                                           const bool &fieldsymmetric)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath << "/components";
  if (Spud::have_option(buffer.str()))
  {                                            
    components = new std::vector<int>;
    serr = Spud::get_option(buffer.str(), *components);            // get the components
    spud_err(buffer.str(), serr);

    if(fieldshape.size()==0)
    {
      tf_err("Requested components of a scalar field.", "Optionpath: %s", buffer.str().c_str());
    }
    else if(fieldshape.size()==1)
    {
      std::vector<int>::iterator max_comp_it =                  
           std::max_element((*components).begin(), (*components).end()); // check the maximum requested component exists
      
      assert(*max_comp_it < fieldshape[0]);
      assert((*components).size() <= fieldshape[0]);
    }
    else if(fieldshape.size()==2)
    {
      std::vector<int>::iterator max_comp_it = 
           std::max_element((*components).begin(), (*components).end()); // check the maximum requested component exists

      if (fieldsymmetric)
      {
        assert(*max_comp_it < fieldshape[0]*fieldshape[1] - ((fieldshape[0]-1)*fieldshape[0])/2);
        assert((*components).size() < fieldshape[0]*fieldshape[1] - ((fieldshape[0]-1)*fieldshape[0])/2);
      }
      else
      {
        assert(*max_comp_it < fieldshape[0]*fieldshape[1]);
        assert((*components).size() <= fieldshape[0]*fieldshape[1]);
      }
    }
    else
    {
      tf_err("Unknown rank.", "fieldshape.size() = %d", fieldshape.size());
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

  buffer.str(""); buffer << optionpath << "/python";
  if (Spud::have_option(buffer.str()))
  {
    std::string pyfunction;
    serr = Spud::get_option(buffer.str(), pyfunction);
    spud_err(buffer.str(), serr);

    if (fieldshape.size()==0)
    {
      value_exp = new PythonExpression(pyfunction);
    }
    else if (fieldshape.size()==1)
    {
      int size = fieldshape[0];
      if (components)
      {
        size = (*components).size();
      }
      value_exp = new PythonExpression(size, pyfunction);
    }
    else if (fieldshape.size()==2)
    {
      if (components)
      {
        value_exp = new PythonExpression((*components).size(), pyfunction);
      }
      else
      {
        value_exp = new PythonExpression(fieldshape, pyfunction);
      }
    }
    else
    {
      tf_err("Unknown rank.", "fieldshape.size() = %d", fieldshape.size());
    }
  }

  buffer.str(""); buffer << optionpath << "/constant";
  if (Spud::have_option(buffer.str()))
  {
    value_const = new double;
    serr = Spud::get_option(buffer.str(), *value_const);
    spud_err(buffer.str(), serr);
  }


}

//*******************************************************************|************************************************************//
// destroy any restrictions we place on the field IS (i.e. region ids, components or boundary_ids)
//*******************************************************************|************************************************************//
void SpudSolverBucket::destroy_field_restrictions_(std::vector<int>* &components,
                                                   std::vector<int>* &region_ids,
                                                   std::vector<int>* &boundary_ids,
                                                   dolfin::Expression* &value_exp,
                                                   double* &value_const)
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
  if (value_exp)
  {
    delete value_exp;
    value_exp = NULL;
  }
  if (value_const)
  {
    delete value_const;
    value_const = NULL;
  }

}

