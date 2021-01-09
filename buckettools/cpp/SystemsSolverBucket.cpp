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


#include "SystemsSolverBucket.h"
#include "SignalHandler.h"
#include "BoostTypes.h"
#include "Bucket.h"
#include "Logger.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SystemsSolverBucket::SystemsSolverBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SystemsSolverBucket::SystemsSolverBucket(const int &solve_location, Bucket* bucket, 
                                         SystemsSolverBucket* systemssolver) : 
                                         solve_location_(solve_location), 
                                         bucket_(bucket), 
                                         systemssolver_(systemssolver)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SystemsSolverBucket::~SystemsSolverBucket()
{
  if(rtol_)
  {
    delete rtol_;
    rtol_ = NULL;
  }
}

//*******************************************************************|************************************************************//
// default filling of the class (i.e. when no user options are specified)
//*******************************************************************|************************************************************//
void SystemsSolverBucket::fill()
{

  fill_base_();

  fill_solvers_();

  fill_solverlists_();

}

//*******************************************************************|************************************************************//
// fill the base data using default values
//*******************************************************************|************************************************************//
void SystemsSolverBucket::fill_base_()
{
  switch(solve_location_)
  {
    case(SOLVE_TIMELOOP):
      name_ = "in_timeloop";
      break;
    case(SOLVE_START):
      name_ =  "at_start";
      break;
    case(SOLVE_DIAGNOSTICS):
      name_ =  "with_diagnostics";
      break;
    default:
      tf_err("Unknown solve location index in SystemsSolverBucket.", "location = %d", solve_location_);
  }

  type_ = "nonlinear_systems_solver";

  rtol_ = NULL;

  atol_ = 1.e-50;

  maxits_ = 1;
 
  minits_ = 0;

  ignore_failures_ = false;

  iteration_count_.reset( new int );
  *iteration_count_ = 0;

  solved_.reset( new bool(false) );                                  // assume the solver hasn't been solved yet

  coupled_ = false;

}

//*******************************************************************|************************************************************//
// fill the list of solvers by iterating over all systems and solvers looking at their solve location flags
//*******************************************************************|************************************************************//
void SystemsSolverBucket::fill_solvers_()
{

  for (SystemBucket_const_it s_it = (*bucket()).systems_begin();
         s_it != (*bucket()).systems_end(); s_it++)
  {
    for (SolverBucket_const_it r_it = (*(*s_it).second).solvers_begin();
           r_it != (*(*s_it).second).solvers_end(); r_it++)
    {
      if ((*(*r_it).second).solve_location() == solve_location())
      {
        register_solver((*r_it).second,                              // make up a, hopefully, unique name
                        (*(*s_it).second).name()+"::"+(*(*r_it).second).name());
        (*(*r_it).second).register_systemssolver(this, unique_name());// reverse registration as well but with our unique-ified name
      }
    }
  }

}

//*******************************************************************|************************************************************//
// fill a list of standard solvers that should be used to calculate the residual of this systems solver
//*******************************************************************|************************************************************//
void SystemsSolverBucket::fill_solverlists_()
{
  
   std::vector<std::string> systemnames;
   GenericSolverBucket_it s_it = solvers_end();
   while (s_it != solvers_begin())                                   // reverse iterate so that solvers that are solved last
   {                                                                 // are registered first and other solvers in the same system
     s_it--;                                                         // aren't used at all
     SystemsSolverBucket_ptr ss_ptr = std::dynamic_pointer_cast<SystemsSolverBucket>((*s_it).second);
     if (ss_ptr)
     {
       const std::vector<SolverBucket_ptr>& tmpresidualsolvers = (*ss_ptr).residualsolvers();
       std::vector<SolverBucket_ptr>::const_iterator ts_it;
       for (ts_it = tmpresidualsolvers.begin(); ts_it != tmpresidualsolvers.end(); ts_it++)
       {                                                             // check we don't have a solver from this system already
         if (std::find(systemnames.begin(), systemnames.end(), (*(**ts_it).system()).name()) == systemnames.end())
         {
           residualsolvers_.push_back(*ts_it);
           systemnames.push_back((*(**ts_it).system()).name());
         }
       }
     }
     else                                                            // treat sub (standard) solvers a bit differently
     {                                                               // as they don't have sub solvers themselves
       SolverBucket_ptr s_ptr = std::dynamic_pointer_cast<SolverBucket>((*s_it).second);
       if (s_ptr)
       {
         if (std::find(systemnames.begin(), systemnames.end(), (*(*s_ptr).system()).name()) == systemnames.end())
         {
           residualsolvers_.push_back(s_ptr);
           systemnames.push_back((*(*s_ptr).system()).name());
         }
       }
       else                                                          // shouldn't hit here (hopefully)
       {
         tf_err("Failed to down cast GenericSolverBucket.", "SystemsSolverBucket name: %s, GenericSolverBucket name: %s", 
                                                                        name().c_str(), (*(*s_it).second).name().c_str());
       }
     }
   }
 
}

//*******************************************************************|************************************************************//
// solve the bilinear system described by the forms in the solver bucket
//*******************************************************************|************************************************************//
bool SystemsSolverBucket::solve()
{
  bool any_solved = false;

  double aerror0 = 0.0;
  *iteration_count_ = 0;

  if (rtol_)
  {
    aerror0 = residual_norm();
    log(INFO, "Entering nonlinear systems iteration: %s", name().c_str());
  }

  while (!complete_iterating_(aerror0))
  {
    (*iteration_count_)++;                                           // increment iteration counter

    for (GenericSolverBucket_const_it s_it = solvers_begin();
                                s_it != solvers_end(); s_it++)
    {
      (*(*s_it).second).set_current_systemssolver(unique_name());
      bool solved = (*(*s_it).second).solve();
      any_solved = any_solved || solved;
      (*(*s_it).second).reset_current_systemssolver();
    }

    //update_nonlinear();
  }

  if (solved_)
  {
    *solved_ = true;
  }

  return any_solved;
  
}

//*******************************************************************|************************************************************//
// a special interface to solve for diagnostics
//*******************************************************************|************************************************************//
bool SystemsSolverBucket::solve_diagnostics(const bool &write_vis,   // a special interface to solve for diagnostics
                                            const bool &write_stat, 
                                            const bool &write_steady, 
                                            const bool &write_det)
{
  bool any_solved = false;

  if (coupled())
  {
    if (!solved())                                                   // if coupled we want to check if any are included in
    {                                                                // diagnostics then solve the whole systems solver
      if( (write_vis    && include_in_visualization()) ||            // but for diagnostics we don't do anything if it's already
          (write_stat   && include_in_statistics())    ||            // been solved for
          (write_steady && include_in_steadystate())   ||
          (write_det    && include_in_detectors())        )
      {
        any_solved = solve();
      }
    }
  }
  else                                                               // if not coupled then we just solve those systems solvers
  {                                                                  // that are required for the diagnostics
    for (GenericSolverBucket_const_it s_it = solvers_begin();
                                s_it != solvers_end(); s_it++)
    {
      if (!(*(*s_it).second).solved())
      {
        if( (write_vis    && (*(*s_it).second).include_in_visualization()) ||
            (write_stat   && (*(*s_it).second).include_in_statistics())    ||
            (write_steady && (*(*s_it).second).include_in_steadystate())   ||
            (write_det    && (*(*s_it).second).include_in_detectors())        )
        {
          bool solved = (*(*s_it).second).solve();
          any_solved = any_solved || solved;
        }
      }
    }
  }

  return any_solved;
}

//*******************************************************************|************************************************************//
// return the l2 norm of the residual of all the systems being solved
//*******************************************************************|************************************************************//
double SystemsSolverBucket::residual_norm()
{
  double norm = 0.0;

  std::vector<SolverBucket_ptr>::const_iterator s_it;
  for (s_it = residualsolvers_.begin(); s_it != residualsolvers_.end(); s_it++)
  {
    norm += std::pow((**s_it).residual_norm(), 2.0);
  }

  norm = std::sqrt(norm);

  return norm;
}

//*******************************************************************|************************************************************//
// update the solver at the end of a timestep
//*******************************************************************|************************************************************//
void SystemsSolverBucket::resetcalculated()
{
  if (solved_)
  {
    *solved_ = false;
  }
  *iteration_count_ = 0;                                             // an iteration counter

  for (GenericSolverBucket_it s_it = solvers_begin(); 
                              s_it != solvers_end(); s_it++)
  {
    (*(*s_it).second).resetcalculated();
  }
}

//*******************************************************************|************************************************************//
// return the number of nonlinear iterations taken
//*******************************************************************|************************************************************//
const int SystemsSolverBucket::iteration_count() const
{
  return *iteration_count_;
}

//*******************************************************************|************************************************************//
// set the number of nonlinear iterations taken
//*******************************************************************|************************************************************//
void SystemsSolverBucket::iteration_count(const int &it)
{
  *iteration_count_ = it;
}

//*******************************************************************|************************************************************//
// return a string that concatenates the names of parent systems solver buckets producing a unique name
//*******************************************************************|************************************************************//
std::string SystemsSolverBucket::unique_name() const
{
  std::string parent_name="";
  if (systemssolver())
  {
    parent_name = (*systemssolver()).unique_name()+"_";
  }
  return parent_name+name();
}

//*******************************************************************|************************************************************//
// a utility function to set the nonlinear systems iterations in output file basenames
//*******************************************************************|************************************************************//
std::string SystemsSolverBucket::iterations_str() const
{
  std::stringstream buffer;
  buffer.str("");
  if (systemssolver())
  {
    buffer << (*systemssolver()).iterations_str();
  }
  if (iterative())
  {
    buffer << "_" << iteration_count();
  }
  return buffer.str();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemsSolverBucket::register_solver(GenericSolverBucket_ptr solver, const std::string &name)
{
  // First check if a solver with this name already exists
  GenericSolverBucket_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it != solvers_.get<om_key_hash>().end())
  {
    tf_err("GenericSolverBucket already exists in systems solver.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    solvers_.insert(om_item<const std::string,GenericSolverBucket_ptr>(name, solver)); // if not then insert it into the solvers_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
GenericSolverBucket_ptr SystemsSolverBucket::fetch_solver(const std::string &name)
{
  GenericSolverBucket_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it == solvers_.get<om_key_hash>().end())
  {
    tf_err("GenericSolverBucket does not exist in systems solver.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does return a (boost shared) pointer to the solver
  }
}

//*******************************************************************|************************************************************//
// return a const (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
const GenericSolverBucket_ptr SystemsSolverBucket::fetch_solver(const std::string &name) const
{
  GenericSolverBucket_const_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it == solvers_.get<om_key_hash>().end())
  {
    tf_err("GenericSolverBucket does not exist in systems solvers solver.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does return a (boost shared) pointer to the solver
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the solvers_ map
//*******************************************************************|************************************************************//
GenericSolverBucket_it SystemsSolverBucket::solvers_begin()
{
  return solvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the solvers_ map
//*******************************************************************|************************************************************//
GenericSolverBucket_const_it SystemsSolverBucket::solvers_begin() const
{
  return solvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
GenericSolverBucket_it SystemsSolverBucket::solvers_end()
{
  return solvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
GenericSolverBucket_const_it SystemsSolverBucket::solvers_end() const
{
  return solvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the systems solver bucket
//*******************************************************************|************************************************************//
const std::string SystemsSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemsSolverBucket " << name() << std::endl;
  indent++;
  s << solvers_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the forms in the solver bucket
//*******************************************************************|************************************************************//
const std::string SystemsSolverBucket::solvers_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( GenericSolverBucket_const_it s_it = solvers_begin(); 
                            s_it != solvers_end(); s_it++ )
  {
    s << indentation << "GenericSolver " << (*s_it).first << std::endl;
    s << (*(*s_it).second).str(indent+1);
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// checkpoint the solverbucket
//*******************************************************************|************************************************************//
void SystemsSolverBucket::checkpoint()
{
  checkpoint_options_();
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this systems solver has fields to be included in visualization output
//*******************************************************************|************************************************************//
const bool SystemsSolverBucket::include_in_visualization() const
{
  bool include = false;
  std::vector<SolverBucket_ptr>::const_iterator sol_it;
  for (sol_it = residualsolvers().begin(); sol_it != residualsolvers().end(); sol_it++)
  {
    const SystemBucket* p_sys = (**sol_it).system();
    include = (*p_sys).include_in_visualization();
    if (include)
    {
      break;
    }
  }
  return include;
}
    
//*******************************************************************|************************************************************//
// return a boolean indicating if this systems solver has fields to be included in statistics output
//*******************************************************************|************************************************************//
const bool SystemsSolverBucket::include_in_statistics() const
{
  bool include = false;
  std::vector<SolverBucket_ptr>::const_iterator sol_it;
  for (sol_it = residualsolvers().begin(); sol_it != residualsolvers().end(); sol_it++)
  {
    const SystemBucket* p_sys = (**sol_it).system();
    include = (*p_sys).include_in_statistics();
    if (include)
    {
      break;
    }
  }
  return include;
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this systems solver has fields to be included in steady state output
//*******************************************************************|************************************************************//
const bool SystemsSolverBucket::include_in_steadystate() const
{
  bool include = false;
  std::vector<SolverBucket_ptr>::const_iterator sol_it;
  for (sol_it = residualsolvers().begin(); sol_it != residualsolvers().end(); sol_it++)
  {
    const SystemBucket* p_sys = (**sol_it).system();
    include = (*p_sys).include_in_steadystate();
    if (include)
    {
      break;
    }
  }
  return include;
}
   
//*******************************************************************|************************************************************//
// return a boolean indicating if this systems solver has fields to be included in detectors output
//*******************************************************************|************************************************************//
const bool SystemsSolverBucket::include_in_detectors() const
{
  bool include = false;
  std::vector<SolverBucket_ptr>::const_iterator sol_it;
  for (sol_it = residualsolvers().begin(); sol_it != residualsolvers().end(); sol_it++)
  {
    const SystemBucket* p_sys = (**sol_it).system();
    include = (*p_sys).include_in_detectors();
    if (include)
    {
      break;
    }
  }
  return include;
}


//*******************************************************************|************************************************************//
// return a boolean indicating if the solver has finished iterating or not
//*******************************************************************|************************************************************//
bool SystemsSolverBucket::complete_iterating_(const double &aerror0)
{
  bool completed = true;
  
  if (rtol_)
  {
    double aerror = residual_norm();
    double rerror;
    if (aerror0 == 0.0)
    {
      rerror = aerror;
    }
    else
    {
      rerror = aerror/aerror0;
    }

    log(INFO, "  %u Nonlinear Systems (%s) Residual Norm (absolute, relative) = %g, %g\n", 
                                    iteration_count(), name().c_str(), aerror, rerror);

    if(convfile_)
    {
      (*convfile_).write_data(aerror);
    }

    if(write_convvis_)
    {
      std::stringstream buffer;
      buffer.str(""); buffer << (*bucket()).output_basename() << "_" << unique_name() << "_" << (*bucket()).timestep_count();
      if (systemssolver())
      {
        buffer << (*systemssolver()).iterations_str();
      }
      buffer << "_systemssolver.xdmf";
      XDMFFile_ptr convvis_file( new dolfin::XDMFFile((*(*(**(residualsolvers_.begin())).system()).mesh()).mpi_comm(), buffer.str()) );
      // making an assumption here that all meshes share the same comm!
      bool append = iteration_count()!=0;

      std::vector<SolverBucket_ptr>::const_iterator s_it;
      for (s_it = residualsolvers_.begin(); s_it != residualsolvers_.end(); s_it++)
      {
        // all fields and residuals get output in this debugging output
        // regardless of whether they're included in standard output
        for (FunctionBucket_it f_it = (*(**s_it).system()).fields_begin();
                               f_it != (*(**s_it).system()).fields_end(); f_it++)
        {
          (*(*f_it).second).write_checkpoint(convvis_file, "iterated", (double)iteration_count(),
                                             append, (*(**s_it).system()).name()+"::Iterated"+(*(*f_it).second).name());
          append=true;
          (*(*f_it).second).write_checkpoint(convvis_file, "residual", (double)iteration_count(),
                                             append, (*(**s_it).system()).name()+"::Residual"+(*(*f_it).second).name());
        }

        // including coefficients here is just a niceity but some
        // coefficients aren't suitable for visualization so
        // only output them if we've asked for them in the normal
        // output
        for (FunctionBucket_it c_it = (*(**s_it).system()).coeffs_begin();
                               c_it != (*(**s_it).system()).coeffs_end(); c_it++)
        {
          if ((*(*c_it).second).include_in_visualization())
          {
            (*(*c_it).second).write_checkpoint(convvis_file, "iterated", (double)iteration_count(),
                                               append, (*(**s_it).system()).name()+"::"+(*(*c_it).second).name());
            append=true;
          }
        }

      }
    }

    completed = ((rerror <= *rtol_ || 
                  aerror <= atol_ || 
                  iteration_count() >= maxits_) 
                 && iteration_count() >= minits_);

    if (iteration_count() == maxits_ && rerror > *rtol_ && aerror > atol_)
    {
      log(WARNING, "it = %d, maxits_ = %d", iteration_count(), maxits_);
      log(WARNING, "rerror = %.12e, rtol_ = %.12e", rerror, *rtol_);
      log(WARNING, "aerror = %.12e, atol_ = %.12e", aerror, atol_);
      if (ignore_failures_)
      {
        log(WARNING, "Ignoring: Nonlinear system failure.");
      }
      else
      {
        tf_fail("Nonlinear systems failed to converge.", "Iteration count, relative error or absolute error too high.");
      }
    }
  }
  else
  {
    completed = iteration_count() >= maxits_;
  }

  if (!completed && (*bucket()).walltime_limit_ptr())
  {
    if ((*bucket()).elapsed_walltime() >= (*bucket()).walltime_limit())
    {
      tf_fail("Nonlinear systems failed to converge.", "Walltime limit reached.");
    }
  }

  if ((*(*SignalHandler::instance()).return_handler(SIGINT)).received())
  {
    log(ERROR, "SIGINT received, terminating nonlinear systems iteration.");
    completed = true;
  }

  return completed;

}

    
//*******************************************************************|************************************************************//
// virtual checkpointing of options
//*******************************************************************|************************************************************//
void SystemsSolverBucket::checkpoint_options_()
{
                                                                     // do nothing - normally we would error but we do declare some
                                                                     // none derived SystemsSolverBuckets so those shouldn't do
                                                                     // anything
}
