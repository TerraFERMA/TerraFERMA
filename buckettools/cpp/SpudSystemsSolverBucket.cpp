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


#include "SpudSystemsSolverBucket.h"
#include "Bucket.h"
#include "BoostTypes.h"
#include "SpudBase.h"
#include "Logger.h"
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudSystemsSolverBucket::SpudSystemsSolverBucket(const std::string &optionpath, 
                                                 const int &solve_location,
                                                 Bucket* bucket, 
                                                 SystemsSolverBucket* systemssolver) : 
                                                 optionpath_(optionpath), 
                                                 SystemsSolverBucket(solve_location, bucket, systemssolver)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudSystemsSolverBucket::~SpudSystemsSolverBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// fill the solver bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSystemsSolverBucket::fill()
{

  fill_base_();                                                      // fill the base solver data: type, name etc.

  fill_solvers_();                                                   // fill the map of sub systems solvers and sub (standard) solvers

  fill_solverlists_();                                               // fill a list of solvers to be used in the residual calculation

}

//*******************************************************************|************************************************************//
// initialize any diagnostic output from the systems solver
//*******************************************************************|************************************************************//
void SpudSystemsSolverBucket::initialize_diagnostics()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::string solver_optionpath;

  std::string l_optionpath = optionpath();                           // the highest level solver doesn't
  buffer.str(""); buffer << optionpath() << "/type";                 // have a type member (compensate for that
  if (Spud::have_option(buffer.str()))                               // here)
  {
    l_optionpath = buffer.str();
  }

  buffer.str(""); buffer << l_optionpath << "/monitors/convergence_file";
  if (Spud::have_option(buffer.str()))
  {
    convfile_.reset( new SystemsConvergenceFile((*bucket()).output_basename()+"_"+unique_name()+".conv",
                                  (*(*(*bucket()).meshes_begin()).second).mpi_comm(),
                                  bucket(), this) );
  }

  for (GenericSolverBucket_it s_it = solvers_begin(); 
                              s_it != solvers_end(); s_it++)
  {
    SpudSystemsSolverBucket_ptr sss_ptr = std::dynamic_pointer_cast<SpudSystemsSolverBucket>((*s_it).second);
    if (sss_ptr)
    {
      (*sss_ptr).initialize_diagnostics();
    }
  }
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the systems solver bucket
//*******************************************************************|************************************************************//
const std::string SpudSystemsSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemsSolverBucket " << name() << " (" << 
                                    optionpath() << ")" << std::endl;
  indent++;
  s << solvers_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// fill the solver bucket base data assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSystemsSolverBucket::fill_base_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // solver name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  std::string l_optionpath = optionpath();                           // the highest level solver doesn't
  buffer.str(""); buffer << optionpath() << "/type";                 // have a type member (compensate for that
  if (Spud::have_option(buffer.str()))                               // here)
  {
    l_optionpath = buffer.str();

    buffer.str(""); buffer << l_optionpath << "/name";               // solver type (as string)
    serr = Spud::get_option(buffer.str(), type_);                    // FIXME: add conversion to enum here
    spud_err(buffer.str(), serr);
  }
  else
  {
    type_ = "nonlinear_systems_solver";
  }

  buffer.str(""); buffer << l_optionpath << "/relative_error";       // relative nonlinear error
  if (Spud::have_option(buffer.str()))
  {
    rtol_ = new double;
    serr = Spud::get_option(buffer.str(), *rtol_); 
    spud_err(buffer.str(), serr);

    coupled_ = true;                                                 // iterative solvers are automatically assumed to be coupled
  }
  else
  {
    buffer.str(""); buffer << l_optionpath << "/coupled";
    coupled_ = Spud::have_option(buffer.str());

    rtol_ = NULL;
  }

  buffer.str(""); buffer << l_optionpath << "/absolute_error";       // absolute nonlinear error
  serr = Spud::get_option(buffer.str(), atol_, 1.e-50); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << l_optionpath << "/max_iterations";       // maximum number of nonlinear iterations
  serr = Spud::get_option(buffer.str(), maxits_, 1); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << l_optionpath << "/min_iterations";       // minimum number of nonlinear iterations (only applies to
  serr = Spud::get_option(buffer.str(), minits_, 0);                 // picard solver types)
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << l_optionpath << 
                                "/ignore_all_solver_failures";
  ignore_failures_ = Spud::have_option(buffer.str());

  iteration_count_.reset( new int );
  *iteration_count_ = 0;

  solved_.reset( new bool(false) );                                  // assume the solver hasn't been solved yet

}

//*******************************************************************|************************************************************//
// fill the list of solvers by iterating over the user options (assuming the buckettools schema)
//*******************************************************************|************************************************************//
void SpudSystemsSolverBucket::fill_solvers_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::string solver_optionpath;

  std::string l_optionpath = optionpath();                           // the highest level solver doesn't
  buffer.str(""); buffer << optionpath() << "/type";                 // have a type member (compensate for that
  if (Spud::have_option(buffer.str()))                               // here)
  {
    l_optionpath = buffer.str();
  }

  buffer.str(""); buffer << l_optionpath << "/solver";
  int nsolvers = Spud::option_count(buffer.str());
  for (uint i=0; i<nsolvers; i++)                                    // loop over the solvers listed
  {
    buffer.str(""); buffer << l_optionpath << "/solver[" << i << "]";
    solver_optionpath = buffer.str();
    
    std::string type;
    buffer.str(""); buffer << solver_optionpath << "/type/name";
    serr = Spud::get_option(buffer.str(), type);
    spud_err(buffer.str(), serr);

    if (type=="nonlinear_systems_solver")
    {
      SpudSystemsSolverBucket_ptr solver( new SpudSystemsSolverBucket(solver_optionpath, 
                                                                      solve_location(), 
                                                                      bucket(), this) );
      (*solver).fill();                                              // recursively fill the sub systems solver

      register_solver(solver, (*solver).name());
    }
    else if (type=="nonlinear_solver")                               // this is just a reference to a standard system::solver
    {
      std::string systemname;
      buffer.str(""); buffer << solver_optionpath << "/type/system/name";
      serr = Spud::get_option(buffer.str(), systemname);
      spud_err(buffer.str(), serr);

      std::string solvername;
      buffer.str(""); buffer << solver_optionpath << "/type/nonlinear_solver/name";
      serr = Spud::get_option(buffer.str(), solvername);
      spud_err(buffer.str(), serr);

      std::string lname;
      buffer.str(""); buffer << solver_optionpath << "/name";
      serr = Spud::get_option(buffer.str(), lname);
      spud_err(buffer.str(), serr);

      SolverBucket_ptr solver = (*(*bucket()).fetch_system(systemname)).fetch_solver(solvername);

      register_solver(solver, lname);                                // not necessarily the name of the solver
                                                                     // instead the name it was given by the user in this list
      (*solver).register_systemssolver(this, unique_name());         // reverse registration as well but with our unique-ified name
    }
    else
    {
      tf_err("Unknown solver type in nonlinear systems solver.", "SystemsSolverBucket: %s", name().c_str());
    }

  }

}

//*******************************************************************|************************************************************//
// checkpoint the options file
//*******************************************************************|************************************************************//
void SpudSystemsSolverBucket::checkpoint_options_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  if (solve_location() == SOLVE_START)                               // do not solve again if this system was only meant to be
  {                                                                  // solved once
    std::string l_optionpath = optionpath();
    buffer.str(""); buffer << optionpath() << "/type";
    if (Spud::have_option(buffer.str()))
    {
      l_optionpath = buffer.str();
    }

    buffer.str(""); buffer << l_optionpath << "/solver";
    int nsolvers = Spud::option_count(buffer.str());
    for (uint i=0; i<nsolvers; i++)
    {
      buffer.str(""); buffer << l_optionpath << "/solver[" << i << "]";
      serr = Spud::delete_option(buffer.str());
      spud_err(buffer.str(), serr);
    }
  }
}

