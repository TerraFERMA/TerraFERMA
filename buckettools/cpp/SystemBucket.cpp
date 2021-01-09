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
#include "Bucket.h"
#include "SystemBucket.h"
#include "FunctionBucket.h"
#include "SolverBucket.h"
#include "SystemsSolverBucket.h"
#include "FunctionalBucket.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SystemBucket::SystemBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SystemBucket::SystemBucket(Bucket* bucket) : bucket_(bucket)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SystemBucket::~SystemBucket()
{
}

//*******************************************************************|************************************************************//
// evaluate the fields at their initial values
//*******************************************************************|************************************************************//
void SystemBucket::evaluate_initial_fields()
{
  log(INFO, "Evaluating initial fields for system %s", name().c_str());

  if (fields_size()>0)
  {
    apply_ic_();                                                     // apply the initial condition to the system function
    apply_bcs_();                                                    // apply the boundary conditions we just collected
  }

}

//*******************************************************************|************************************************************//
// loop over the ordered solver buckets in the bucket, calling solve on each of them
//*******************************************************************|************************************************************//
bool SystemBucket::solve(const int &location, const bool force)
{
  std::vector<int> locations(1, location);
  return solve(locations, force);
}

//*******************************************************************|************************************************************//
// loop over the ordered solver buckets in the bucket, calling solve on each of them
//*******************************************************************|************************************************************//
bool SystemBucket::solve(const std::vector<int> &locations, const bool force)
{
  bool any_solved = false;

  // the introduction of systems solvers makes this logic quite complicated...
  for (i_SystemsSolverBucket_it ss_it = (*bucket()).systemssolvers_begin();
                  ss_it != (*bucket()).systemssolvers_end(); ss_it++)
  {                                                                  // loop over the highest level systems solvers in the bucket
    bool solve = true;
    if (!locations.empty())                                          // skip any system solve that isn't at the requested location
    {
      solve = std::find(locations.begin(), locations.end(), (*ss_it).first) != locations.end();
    }
    if (!solve)
    {
      continue;
    }

    const std::vector<SolverBucket_ptr>& residualsolvers = (*(*ss_it).second).residualsolvers();
    bool solver_present = false;                                     // also skip any system solver that doesn't contain a solver 
    for (SolverBucket_const_it s_it = solvers_begin(); s_it != solvers_end(); s_it++) // from this system
    {
      if (std::find(residualsolvers.begin(), residualsolvers.end(), (*s_it).second) != residualsolvers.end())
      {
        solver_present = true;
        break;
      }
    }
    if (!solver_present)
    {
      continue;
    }

    // if we made it this far then we are in a systems solver that meets the location criterion and has something to do with this
    // sytem
    if ((*(*ss_it).second).coupled())                                // this is simple if we're coupled, we just solve the whole
    {                                                                // systems solver
      if (!(*(*ss_it).second).solved() || any_solved || force)       // but only if it hasn't been solved before or a previous
      {                                                              // solver has been solved or we're forcing it
        bool solved = (*(*ss_it).second).solve();
        any_solved = any_solved || solved;
      }
    }
    else                                                             // if we're not coupled then we only want to solve those
    {                                                                // sub generic solvers that are relevant to this sytem
      for (GenericSolverBucket_const_it gs_it = (*(*ss_it).second).solvers_begin();// loop over the sub generic solvers
                                  gs_it != (*(*ss_it).second).solvers_end(); gs_it++)
      {
        SystemsSolverBucket_ptr ss_ptr = std::dynamic_pointer_cast<SystemsSolverBucket>((*gs_it).second);
        if (ss_ptr)                                                  // if this is a systems solver then
        {
          const std::vector<SolverBucket_ptr>& subresidualsolvers = (*ss_ptr).residualsolvers();
          bool subsolver_present = false;                            // test if it touches solvers in this sytem and skip if
          for (SolverBucket_const_it s_it = solvers_begin(); s_it != solvers_end(); s_it++) // it doesn't
          {
            if (std::find(subresidualsolvers.begin(), subresidualsolvers.end(), (*s_it).second) != subresidualsolvers.end())
            {
              subsolver_present = true;
              break;
            }
          }
          if (!subsolver_present)
          {
            continue;
          }

          if (!(*ss_ptr).solved() || any_solved || force)
          {
            bool solved = (*ss_ptr).solve();
            any_solved = any_solved || solved;
          }
        }
        else
        {
          SolverBucket_ptr s_ptr = std::dynamic_pointer_cast<SolverBucket>((*gs_it).second);
          if (s_ptr)                                                 // if this is a standard solver then
          {                                                          // just test that it's from our system
            bool subsolver_present = false;
            for (SolverBucket_const_it s_it = solvers_begin(); s_it != solvers_end(); s_it++)
            {
              if ((*s_it).second == s_ptr)
              {
                subsolver_present = true;
                break;
              }
            }
            if (!subsolver_present)
            {
              continue;
            }

            if (!(*s_ptr).solved() || any_solved || force)
            {
              bool solved = (*s_ptr).solve();
              any_solved = any_solved || solved;
            }
          }
          else
          {
            tf_err("Failed to down cast GenericSolverBucket.", "SystemBucket name: %s, GenericSolverBucket name: %s", 
                                                                           name().c_str(), (*(*gs_it).second).name().c_str());
          }
        }
      }
    }
  }

  return any_solved;
}

//*******************************************************************|************************************************************//
// update the timelevels of the system function
//*******************************************************************|************************************************************//
void SystemBucket::update()
{
  if (function_)
  {
    (*(*oldfunction_).vector()) = (*(*function_).vector());          // update the oldfunction to the new function value
  }
  
                                                                     // fields share a vector with the system function so no need to
                                                                     // update them...
  for (FunctionBucket_it f_it = coeffs_begin();                      // so just loop over coefficients again to update any coefficient
                           f_it != coeffs_end(); f_it++)             // functions, constant functionals or statistic functionals
  {
    (*(*f_it).second).update();
  }

  for (FunctionalBucket_it f_it = functionals_begin(); 
                                  f_it != functionals_end(); f_it++)
  {
    (*(*f_it).second).update();
  }

  resetcalculated();                                                 // NOTE: this must be done LAST as some functionals may
                                                                     // actually be calculated by this routine to ensure any steady 
                                                                     // state check that uses them is accurate
}

//*******************************************************************|************************************************************//
// update the potentially time dependent components of the system
//*******************************************************************|************************************************************//
void SystemBucket::update_timedependent()
{

  for (FunctionBucket_it f_it = coeffs_begin();                      // loop over coefficients again to update any function coefficients
                           f_it != coeffs_end(); f_it++)
  {
    (*(*f_it).second).update_timedependent();                        // this does nothing to non function coefficients
  }

}

//*******************************************************************|************************************************************//
// update the potentially nonlinear components of the system
//*******************************************************************|************************************************************//
void SystemBucket::update_nonlinear()
{

  for (FunctionBucket_it f_it = coeffs_begin();                      // loop over coefficients again to update any constant
                           f_it != coeffs_end(); f_it++)             // functionals
  {
    (*(*f_it).second).update_nonlinear();                            // this does nothing to non constant functionals
  }

}

//*******************************************************************|************************************************************//
// return the maximum change in the requested fields of the system
//*******************************************************************|************************************************************//
const double SystemBucket::maxchange()
{
  double maxchange = 0.0;

  updatechange();

  for (FunctionBucket_it f_it = fields_begin(); 
                                  f_it != fields_end(); f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())
    {
      double fieldchange = (*(*f_it).second).change();
      log(DBG, "    steady state fieldchange = %f", fieldchange);
      maxchange = std::max( fieldchange, maxchange );
    }
  }

  for (FunctionalBucket_it f_it = functionals_begin(); 
                                  f_it != functionals_end(); f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())
    {
      double functionalchange = (*(*f_it).second).change();
      log(DBG, "      steady state functionalchange = %f", functionalchange);
      maxchange = std::max( functionalchange, maxchange );
    }
  }
  return maxchange;
}

//*******************************************************************|************************************************************//
// update the change function
//*******************************************************************|************************************************************//
void SystemBucket::updatechange()
{
  if (!*change_calculated_ && changefunction_)                       // changefunction_ won't be associated for systems with no fields
  {
    (*(*changefunction_).vector()) = (*(*function_).vector());       // before updating the oldfunction to the new values
    (*(*changefunction_).vector()) -= (*(*oldfunction_).vector());   // update the change in the fields over this timesteps

    *change_calculated_ = true;
  }
}

//*******************************************************************|************************************************************//
// reset the calculated booleans
//*******************************************************************|************************************************************//
void SystemBucket::resetcalculated()
{
  if (change_calculated_)
  {
    *change_calculated_ = false;
  }

  for (FunctionBucket_it f_it = fields_begin();
                           f_it != fields_end(); f_it++)
  {
    (*(*f_it).second).resetcalculated();
  }

  for (FunctionBucket_it f_it = coeffs_begin();
                           f_it != coeffs_end(); f_it++)
  {
    (*(*f_it).second).resetcalculated();
  }

  for (SolverBucket_it s_it = solvers_begin();
                         s_it != solvers_end(); s_it++)
  {
    (*(*s_it).second).resetcalculated();                                      // reset the solved_ indicator to false for the next timestep
  }

  for (FunctionalBucket_it f_it = functionals_begin(); 
                                  f_it != functionals_end(); f_it++)
  {
    (*(*f_it).second).resetcalculated();
  }

}

//*******************************************************************|************************************************************//
// cap the value of the fields in this system
//*******************************************************************|************************************************************//
void SystemBucket::postprocess_values()
{
  for (FunctionBucket_it f_it = fields_begin();
                                  f_it != fields_end(); f_it++)
  {
    (*(*f_it).second).postprocess_values();
  }
}

//*******************************************************************|************************************************************//
// return the l2 norm of the residual of the last solver that meets the location requirement
//*******************************************************************|************************************************************//
double SystemBucket::residual_norm(const int &location)
{
  std::vector<int> locations(1, location);
  return residual_norm(locations);
}

//*******************************************************************|************************************************************//
// return the l2 norm of the residual of the last solver that meets the location requirement
//*******************************************************************|************************************************************//
double SystemBucket::residual_norm(const std::vector<int> &locations)
{
  double norm = 0.0;

  if (!solvers_.empty())
  {
    bool calculate_norm = true;

    SolverBucket_it s_it = solvers_end();
    while (s_it != solvers_begin())
    {
      s_it--;
      if (!locations.empty())
      {
        calculate_norm = std::find(locations.begin(), locations.end(), (*(*s_it).second).solve_location()) != locations.end();
      }
      if (calculate_norm)
      {
        break;
      }
    }

    if (calculate_norm)
    {
      norm = (*(*s_it).second).residual_norm();
    }
  }

  return norm;
}

//*******************************************************************|************************************************************//
// return a std::vector listing the solve locations of this system's solvers
//*******************************************************************|************************************************************//
const std::vector<int> SystemBucket::solve_locations() const
{
  std::vector<int> locations;
  for (SolverBucket_it s_it = solvers_begin();
                         s_it != solvers_end();
                         s_it++)
  {
    locations.push_back((*(*s_it).second).solve_location());
  }

  return locations;
}

//*******************************************************************|************************************************************//
// attach all coefficients to the functionals and solver forms
//*******************************************************************|************************************************************//
void SystemBucket::initialize_forms()
{
  log(INFO, "Attaching coeffs for system %s", name().c_str());
  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end(); 
                                                              f_it++)// attach functions to coefficients constant functionals
  {
    (*(*f_it).second).attach_form_coeffs();
                                          
  }
  for (SolverBucket_it s_it = solvers_begin(); s_it != solvers_end(); 
                                                              s_it++)// attach functions to the solver forms
  {
    (*(*s_it).second).attach_form_coeffs();
  }
  for (FunctionalBucket_it f_it = functionals_begin(); 
                                   f_it != functionals_end(); f_it++)// attach functions to the functional forms
  {
    (*(*f_it).second).attach_form_coeffs();
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a generic function, which one depends on the time
//*******************************************************************|************************************************************//
const Function_ptr SystemBucket::function_ptr(const double_ptr time) const
{
  if (time == (*bucket()).current_time_ptr())
  {
    return iteratedfunction_;
  }
  else if (time == (*bucket()).old_time_ptr())
  {
    return oldfunction_;
  }
  else
  {
    tf_err("Unknown time pointer when returning function in function_ptr.", "SystemBucket: %s", name_.c_str());
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_field(FunctionBucket_ptr field, const std::string &name)
{
  FunctionBucket_hash_it f_it = fields_.get<om_key_hash>().find(name);                       // check if name already exists
  if (f_it != fields_.get<om_key_hash>().end())
  {
    tf_err("Field already exists in system.", "Field name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    fields_.insert(om_item<const std::string, FunctionBucket_ptr>(name,field));                                           // if not, add it to the fields_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_field(const std::string &name)
{
  FunctionBucket_hash_it f_it = fields_.get<om_key_hash>().find(name);                       // check if name already exists
  if (f_it == fields_.get<om_key_hash>().end())
  {
    tf_err("Field does not exist in system.", "Field name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does return a (boost shared) pointer to the field
  }
}

//*******************************************************************|************************************************************//
// return a constant (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
const FunctionBucket_ptr SystemBucket::fetch_field(const std::string &name) const
{
  FunctionBucket_hash_it f_it = fields_.get<om_key_hash>().find(name);                       // check if name already exists
  if (f_it == fields_.get<om_key_hash>().end())
  {
    tf_err("Field does not exist in system.", "Field name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does return a constant (boost shared) pointer to the
                                                                     // field
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::fields_begin()
{
  return fields_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::fields_begin() const
{
  return fields_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::fields_end()
{
  return fields_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::fields_end() const
{
  return fields_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return the number of fields registered
//*******************************************************************|************************************************************//
const int SystemBucket::fields_size() const
{
  return fields_.size();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_coeff(FunctionBucket_ptr coeff, const std::string &name)
{
  FunctionBucket_hash_it f_it = coeffs_.get<om_key_hash>().find(name);                       // check if name already exists
  if (f_it != coeffs_.get<om_key_hash>().end())
  {
    tf_err("Coefficient already exists in system.", "Coefficient name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    coeffs_.insert(om_item<const std::string,FunctionBucket_ptr>(name,coeff));          // if it doesn't, register the function bucket
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_coeff(const std::string &name)
{
  FunctionBucket_hash_it f_it = coeffs_.get<om_key_hash>().find(name);                       // check if the name already exists
  if (f_it == coeffs_.get<om_key_hash>().end())
  {
    tf_err("Coefficient does not exist in system.", "Coefficient name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return the coefficient
  }
}

//*******************************************************************|************************************************************//
// return a constant (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
const FunctionBucket_ptr SystemBucket::fetch_coeff(const std::string &name) const
{
  FunctionBucket_const_hash_it f_it = coeffs_.get<om_key_hash>().find(name);                 // check if the name already exists
  if (f_it == coeffs_.get<om_key_hash>().end())
  {
    tf_err("Coefficient does not exist in system.", "Coefficient name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return the coefficient
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::coeffs_begin()
{
  return coeffs_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::coeffs_begin() const
{
  return coeffs_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::coeffs_end()
{
  return coeffs_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::coeffs_end() const
{
  return coeffs_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_solver(SolverBucket_ptr solver, const std::string &name)
{
  // First check if a solver with this name already exists
  SolverBucket_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it != solvers_.get<om_key_hash>().end())
  {
    tf_err("SolverBucket already exists in system.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    solvers_.insert(om_item<const std::string,SolverBucket_ptr>(name, solver)); // if not then insert it into the solvers_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
SolverBucket_ptr SystemBucket::fetch_solver(const std::string &name)
{
  SolverBucket_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it == solvers_.get<om_key_hash>().end())
  {
    tf_err("SolverBucket does not exist in system.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does return a (boost shared) pointer to the solver
  }
}

//*******************************************************************|************************************************************//
// return a const (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
const SolverBucket_ptr SystemBucket::fetch_solver(const std::string &name) const
{
  SolverBucket_hash_it s_it = solvers_.get<om_key_hash>().find(name);     // check if a solver with this name already exists
  if (s_it == solvers_.get<om_key_hash>().end())
  {
    tf_err("SolverBucket does not exist in system.", "SolverBucket name: %s, System name: %s", name.c_str(), name_.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does return a (boost shared) pointer to the solver
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_it SystemBucket::solvers_begin()
{
  return solvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_const_it SystemBucket::solvers_begin() const
{
  return solvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_it SystemBucket::solvers_end()
{
  return solvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_const_it SystemBucket::solvers_end() const
{
  return solvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a functional form in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_functional(FunctionalBucket_ptr functional, const std::string &name)
{
  FunctionalBucket_hash_it f_it = functionals_.get<om_key_hash>().find(name);                            // check if the name already exists
  if (f_it != functionals_.get<om_key_hash>().end())
  {
    tf_err("Functional already exists in system.", "Functional name: %s, System name: %s", 
           name.c_str(), name_.c_str());
  }
  else
  {
    functionals_.insert(om_item<const std::string, FunctionalBucket_ptr>(name, functional));                                 // if not, insert it into the functionals_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a functional form from the system bucket data maps
//*******************************************************************|************************************************************//
FunctionalBucket_ptr SystemBucket::fetch_functional(const std::string &name)
{
  FunctionalBucket_hash_it f_it = functionals_.get<om_key_hash>().find(name);                            // check if the name already exists
  if (f_it == functionals_.get<om_key_hash>().end())
  {
    tf_err("Functional does not exist in system.", "Functional name: %s, System name: %s", 
           name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a constant (boost shared) pointer to a functional form from the system bucket data maps
//*******************************************************************|************************************************************//
const FunctionalBucket_ptr SystemBucket::fetch_functional(const std::string &name) const
{
  FunctionalBucket_hash_it f_it = functionals_.get<om_key_hash>().find(name);                            // check if the name already exists
  if (f_it == functionals_.get<om_key_hash>().end())
  {
    tf_err("Functional does not exist in system.", "Functional name: %s, System name: %s", 
           name.c_str(), name_.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the functionals_ map
//*******************************************************************|************************************************************//
FunctionalBucket_it SystemBucket::functionals_begin()
{
  return functionals_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the functionals_ map
//*******************************************************************|************************************************************//
FunctionalBucket_const_it SystemBucket::functionals_begin() const
{
  return functionals_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the functionals_ map
//*******************************************************************|************************************************************//
FunctionalBucket_it SystemBucket::functionals_end()
{
  return functionals_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the functionals_ map
//*******************************************************************|************************************************************//
FunctionalBucket_const_it SystemBucket::functionals_end() const
{
  return functionals_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector< std::shared_ptr<const dolfin::DirichletBC> >::iterator SystemBucket::bcs_begin()
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator SystemBucket::bcs_begin() const
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector< std::shared_ptr<const dolfin::DirichletBC> >::iterator SystemBucket::bcs_end()
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator SystemBucket::bcs_end() const
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this system has fields to be included in visualization output
//*******************************************************************|************************************************************//
const bool SystemBucket::include_in_visualization() const
{
  bool include = false;
  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    include = (*(*f_it).second).include_in_visualization() || (*(*f_it).second).include_residual_in_visualization();
    if (include)
    {
      break;
    }
  }

  if (!include)
  {
    for (FunctionBucket_const_it c_it = coeffs_begin(); c_it != coeffs_end(); 
                                                                c_it++)
    {
      include = (*(*c_it).second).include_in_visualization();
      if (include)
      {
        break;
      }
    }
  }

  return include;
}
    
//*******************************************************************|************************************************************//
// return a boolean indicating if this system has fields to be included in statistics output
//*******************************************************************|************************************************************//
const bool SystemBucket::include_in_statistics() const
{
  bool include = false;
  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    include = (*(*f_it).second).include_in_statistics();
    if (include)
    {
      break;
    }
  }

  if (!include)
  {
    for (FunctionalBucket_const_it f_it = functionals_begin(); f_it != functionals_end(); 
                                                                f_it++)
    {
      include = (*(*f_it).second).include_in_statistics();
      if (include)
      {
        break;
      }
    }
  }
  return include;
}
    
//*******************************************************************|************************************************************//
// return a boolean indicating if this system has fields to be included in steady state output
//*******************************************************************|************************************************************//
const bool SystemBucket::include_in_steadystate() const
{
  bool include = false;
  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    include = (*(*f_it).second).include_in_steadystate();
    if (include)
    {
      break;
    }
  }

  if (!include)
  {
    for (FunctionalBucket_const_it f_it = functionals_begin(); f_it != functionals_end(); 
                                                                f_it++)
    {
      include = (*(*f_it).second).include_in_steadystate();
      if (include)
      {
        break;
      }
    }
  }
  return include;
}
   
//*******************************************************************|************************************************************//
// return a boolean indicating if this system has fields to be included in detectors output
//*******************************************************************|************************************************************//
const bool SystemBucket::include_in_detectors() const
{
  bool include = false;
  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    include = (*(*f_it).second).include_in_detectors();
    if (include)
    {
      break;
    }
  }
  return include;
}
    
//*******************************************************************|************************************************************//
// return a string describing the contents of the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << std::endl;
  indent++;
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  s << functionals_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the fields in the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::fields_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionBucket_const_it f_it = fields_begin();              // loop over the fields
                                    f_it != fields_end(); f_it++ )
  {
    s << (*(*f_it).second).str(indent);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the coefficients in the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::coeffs_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionBucket_const_it f_it = coeffs_begin();              // loop over the coefficients
                                  f_it != coeffs_end(); f_it++ )
  {
    s << (*(*f_it).second).str(indent);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the solvers in the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::solvers_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( SolverBucket_const_it s_it = solvers_begin();               // loop over the solvers
                                s_it != solvers_end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
    s << (*(*s_it).second).forms_str(indent+1);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the functionals in the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::functionals_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionalBucket_const_it f_it = functionals_begin(); f_it != functionals_end(); f_it++ )
  {
    s << (*(*f_it).second).str(indent);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// output the system
//*******************************************************************|************************************************************//
void SystemBucket::output()
{

  for (FunctionBucket_it f_it = fields_begin();
                         f_it != fields_end(); f_it++)
  {
    (*(*f_it).second).output();
  }

  for (FunctionBucket_it c_it = coeffs_begin();
                         c_it != coeffs_end(); c_it++)
  {
    (*(*c_it).second).output();
  }

}

//*******************************************************************|************************************************************//
// checkpoint the system
//*******************************************************************|************************************************************//
void SystemBucket::checkpoint(const double_ptr time)
{

  for (FunctionBucket_it f_it = fields_begin();                      // if there's no function then there should be no fields
                         f_it != fields_end(); f_it++)               // so this is a bit redundant outside the above if statement
  {
    (*(*f_it).second).checkpoint(time);
  }

  for (SolverBucket_it s_it = solvers_begin();
                       s_it != solvers_end(); s_it++)
  {
    (*(*s_it).second).checkpoint();
  }

}

//*******************************************************************|************************************************************//
// initialize the system with a combined initial condition (calls eval)
//*******************************************************************|************************************************************//
void SystemBucket::apply_ic_()
{
  for (FunctionBucket_it f_it = fields_begin();                      // if there's no function then there should be no fields
                         f_it != fields_end(); f_it++)               // so this is a bit redundant outside the above if statement
  {
    (*(*f_it).second).apply_ic();
  }
  (*(*iteratedfunction_).vector()) = (*(*oldfunction_).vector());    // set the iterated function vector to the old function vector
  (*(*function_).vector()) = (*(*oldfunction_).vector());            // set the function vector to the old function vector
}

//*******************************************************************|************************************************************//
// apply the vector of system boundary conditions to the system function vectors to ensure consisten initial and boundary conditions
//*******************************************************************|************************************************************//
void SystemBucket::apply_bcs_()
{
  for (std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator     // loop over all the bcs
       b_it = bcs_begin(); b_it != bcs_end(); b_it++)
  {
    (**b_it).apply((*(*oldfunction_).vector()));
    (**b_it).apply((*(*iteratedfunction_).vector()));
    (**b_it).apply((*(*function_).vector()));
  }
}


