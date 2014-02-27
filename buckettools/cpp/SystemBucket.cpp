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


#include "InitialConditionExpression.h"
#include "BoostTypes.h"
#include "SystemBucket.h"
#include "FunctionBucket.h"
#include "SolverBucket.h"
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
  empty_();                                                          // empty the data structures
}

//*******************************************************************|************************************************************//
// loop over the ordered solver buckets in the bucket, calling solve on each of them
//*******************************************************************|************************************************************//
void SystemBucket::solve()
{
  for (int_SolverBucket_const_it s_it = orderedsolvers_begin(); 
                                s_it != orderedsolvers_end(); s_it++)
  {
    (*(*s_it).second).solve();

    cap_values();

    (*(*residualfunction_).vector()) = (*boost::dynamic_pointer_cast< dolfin::GenericVector >((*(*s_it).second).residual_vector()));
    // update_nonlinear...
  }

  if(solved_)
  {
    *solved_ = true;
  }
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
  for (int_FunctionBucket_it f_it = orderedfields_begin();           // except that they contain functionals which need updating
                           f_it != orderedfields_end(); f_it++) 
  {
    (*(*f_it).second).update();
  }

  for (int_FunctionBucket_it f_it = orderedcoeffs_begin();           // also loop over coefficients again to update any coefficient
                           f_it != orderedcoeffs_end(); f_it++)      // functions, constant functionals or statistic functionals
  {
    (*(*f_it).second).update();
  }

  resetcalculated();                                                 // reset the calculated booleans in the system, fields and functionals

  *solved_ = false;                                                  // reset the solved_ indicator to false for the next timestep

}

//*******************************************************************|************************************************************//
// update the potentially time dependent components of the system
//*******************************************************************|************************************************************//
void SystemBucket::update_timedependent()
{

  for (int_FunctionBucket_it f_it = orderedcoeffs_begin();           // loop over coefficients again to update any constant
                           f_it != orderedcoeffs_end(); f_it++)      // functionals
  {
    (*(*f_it).second).update_timedependent();                        // this does nothing to non constant functionals
  }

}

//*******************************************************************|************************************************************//
// update the potentially nonlinear components of the system
//*******************************************************************|************************************************************//
void SystemBucket::update_nonlinear()
{

  for (int_FunctionBucket_it f_it = orderedcoeffs_begin();           // loop over coefficients again to update any constant
                           f_it != orderedcoeffs_end(); f_it++)      // functionals
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
      dolfin::log(dolfin::DBG, "    steady state fieldchange = %f", fieldchange);
      maxchange = std::max( fieldchange, maxchange );
    }

    for (Form_const_it form_it = (*(*f_it).second).functionals_begin();
              form_it != (*(*f_it).second).functionals_end(); form_it++)
    {
      if ((*(*f_it).second).include_functional_in_steadystate((*form_it).first))
      {
        double functionalchange = (*(*f_it).second).functionalchange(form_it);
        dolfin::log(dolfin::DBG, "      steady state functionalchange = %f", functionalchange);
        maxchange = std::max( functionalchange, maxchange );
      }
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
  *change_calculated_ = false;
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
}

//*******************************************************************|************************************************************//
// cap the value of the fields in this system
//*******************************************************************|************************************************************//
void SystemBucket::cap_values()
{
  for (FunctionBucket_it f_it = fields_begin();
                                  f_it != fields_end(); f_it++)
  {
    (*(*f_it).second).cap_values();
  }
}

//*******************************************************************|************************************************************//
// make a partial copy of the provided system bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SystemBucket::copy_diagnostics(SystemBucket_ptr &system, Bucket_ptr &bucket) const
{

  if(!system)
  {
    system.reset( new SystemBucket(&(*bucket)) );
  }

  (*system).name_ = name_;

  (*system).mesh_ = mesh_;

  (*system).functionspace_ = functionspace_;

  (*system).function_ = function_;
  (*system).iteratedfunction_ = iteratedfunction_;
  (*system).oldfunction_ = oldfunction_;

  (*system).changefunction_ = changefunction_;
  (*system).change_calculated_ = change_calculated_;

  (*system).residualfunction_ = residualfunction_;
  (*system).snesupdatefunction_ = snesupdatefunction_;

  (*system).solved_ = solved_;

  for (FunctionBucket_const_it func_it = fields_begin();             // loop over the fields
                           func_it != fields_end(); func_it++)
  {                                                                  
    FunctionBucket_ptr field;                                        // create a new field
    
    (*(*func_it).second).copy_diagnostics(field, system);

    (*system).register_field(field, (*field).name());                // put the field in the bucket
  }                                                                  

  for (FunctionBucket_const_it func_it = coeffs_begin();             // loop over the fields
                           func_it != coeffs_end(); func_it++)
  {                                                                  
    FunctionBucket_ptr coeff;                                        // create a new coeff

    (*(*func_it).second).copy_diagnostics(coeff, system);

    (*system).register_coeff(coeff, (*coeff).name());                // put the coefficient in the bucket
  }                                                                  

  for (SolverBucket_const_it solver_it = solvers_begin();            // loop over the solvers
                      solver_it != solvers_end(); solver_it++)
  {                                                                  
    SolverBucket_ptr solver;                                         // create a new solver
  
    (*(*solver_it).second).copy_diagnostics(solver, system);

    (*system).register_solver(solver, (*solver).name());             // put the coefficient in the bucket
  }                                                                  

}

//*******************************************************************|************************************************************//
// initialize any diagnostic output in this system
//*******************************************************************|************************************************************//
void SystemBucket::initialize_diagnostics() const                    // doesn't allocate anything so can be const
{
  for (SolverBucket_const_it s_it = solvers_begin(); s_it != solvers_end(); s_it++)
  {
    (*(*s_it).second).initialize_diagnostics();
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
    dolfin::error("Unknown time pointer when returning function in genericfunction(time).");
  }
}

//*******************************************************************|************************************************************//
// return the residual associated with the last solver in this system
//*******************************************************************|************************************************************//
const PETScVector_ptr SystemBucket::residual_vector() const
{
  return (*(*(orderedsolvers_.lower_bound((int) solvers_.size()))).second).residual_vector();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_field(FunctionBucket_ptr field, const std::string &name)
{
  FunctionBucket_it f_it = fields_.find(name);                       // check if name already exists
  if (f_it != fields_.end())
  {
    dolfin::error("Field named \"%s\" already exists in system.",    // if it does, issue an error
                                                name.c_str());
  }
  else
  {
    fields_[name] = field;                                           // if not, add it to the fields_ map
    orderedfields_[(int) fields_.size()] = field;                    // and into the orderedfields_ map, assuming that the
                                                                     // insertion order is the order they are to be calculated
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_field(const std::string &name)
{
  FunctionBucket_it f_it = fields_.find(name);                       // check if name already exists
  if (f_it == fields_.end())
  {
    dolfin::error("Field named \"%s\" does not exist in system.",    // if it doesn't, issue an error
                                                    name.c_str());
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
  FunctionBucket_const_it f_it = fields_.find(name);                 // check if name already exists
  if (f_it == fields_.end())
  {
    dolfin::error("Field named \"%s\" does not exist in system.",    // if it doesn't, issue an error
                                                    name.c_str());
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
  return fields_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::fields_begin() const
{
  return fields_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::fields_end()
{
  return fields_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the fields_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::fields_end() const
{
  return fields_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedfields_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_it SystemBucket::orderedfields_begin()
{
  return orderedfields_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedfields_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_const_it SystemBucket::orderedfields_begin() const
{
  return orderedfields_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedfields_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_it SystemBucket::orderedfields_end()
{
  return orderedfields_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedfields_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_const_it SystemBucket::orderedfields_end() const
{
  return orderedfields_.end();
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
  FunctionBucket_it f_it = coeffs_.find(name);                       // check if name already exists
  if (f_it != coeffs_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
              "Coefficient named \"%s\" already exists in system.", 
                                                    name.c_str());
  }
  else
  {
    coeffs_[name] = coeff;                                           // if it doesn't, register the function bucket
    orderedcoeffs_[(int) coeffs_.size()] = coeff;                    // and into the orderedcoeffs_ map, assuming that the
                                                                     // insertion order is the order they are to be calculated
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_coeff(const std::string &name)
{
  FunctionBucket_it f_it = coeffs_.find(name);                       // check if the name already exists
  if (f_it == coeffs_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Coefficient named \"%s\" does not exist in system.", 
                                                    name.c_str());
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
  FunctionBucket_const_it f_it = coeffs_.find(name);                 // check if the name already exists
  if (f_it == coeffs_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Coefficient named \"%s\" does not exist in system.", 
                                                    name.c_str());
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
  return coeffs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::coeffs_begin() const
{
  return coeffs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_it SystemBucket::coeffs_end()
{
  return coeffs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the coeffs_ map
//*******************************************************************|************************************************************//
FunctionBucket_const_it SystemBucket::coeffs_end() const
{
  return coeffs_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedcoeffs_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_it SystemBucket::orderedcoeffs_begin()
{
  return orderedcoeffs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedcoeffs_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_const_it SystemBucket::orderedcoeffs_begin() const
{
  return orderedcoeffs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedcoeffs_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_it SystemBucket::orderedcoeffs_end()
{
  return orderedcoeffs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedcoeffs_ map
//*******************************************************************|************************************************************//
int_FunctionBucket_const_it SystemBucket::orderedcoeffs_end() const
{
  return orderedcoeffs_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_solver(SolverBucket_ptr solver, const std::string &name)
{
  // First check if a solver with this name already exists
  SolverBucket_it s_it = solvers_.find(name);                        // check if this name exists already
  if (s_it != solvers_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
              "SolverBucket named \"%s\" already exists in system.", 
                                                      name.c_str());
  }
  else
  {
    solvers_[name] = solver;                                         // if not then insert it into the solvers_ map
    orderedsolvers_[(int) solvers_.size()] = solver;                 // and into the orderedsolvers_ map, assuming that the
                                                                     // insertion order is the order they are to be solved
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
SolverBucket_ptr SystemBucket::fetch_solver(const std::string &name)
{
  SolverBucket_it s_it = solvers_.find(name);                        // check if name already exists
  if (s_it == solvers_.end())
  {
    dolfin::error("SolverBucket named \"%s\" does not exist in system.",// if it doesn't, issue an error
                                                    name.c_str());
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
  SolverBucket_const_it s_it = solvers_.find(name);                        // check if name already exists
  if (s_it == solvers_.end())
  {
    dolfin::error("SolverBucket named \"%s\" does not exist in system.",// if it doesn't, issue an error
                                                    name.c_str());
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
  return solvers_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_const_it SystemBucket::solvers_begin() const
{
  return solvers_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_it SystemBucket::solvers_end()
{
  return solvers_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the solvers_ map
//*******************************************************************|************************************************************//
SolverBucket_const_it SystemBucket::solvers_end() const
{
  return solvers_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedsolvers_ map
//*******************************************************************|************************************************************//
int_SolverBucket_it SystemBucket::orderedsolvers_begin()
{
  return orderedsolvers_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedsolvers_ map
//*******************************************************************|************************************************************//
int_SolverBucket_const_it SystemBucket::orderedsolvers_begin() const
{
  return orderedsolvers_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedsolvers_ map
//*******************************************************************|************************************************************//
int_SolverBucket_it SystemBucket::orderedsolvers_end()
{
  return orderedsolvers_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedsolvers_ map
//*******************************************************************|************************************************************//
int_SolverBucket_const_it SystemBucket::orderedsolvers_end() const
{
  return orderedsolvers_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the dirichletbcs_ vector
//*******************************************************************|************************************************************//
std::vector< const dolfin::DirichletBC* >::iterator SystemBucket::dirichletbcs_begin()
{
  return dirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the dirichletbcs_ vector
//*******************************************************************|************************************************************//
std::vector< const dolfin::DirichletBC* >::const_iterator SystemBucket::dirichletbcs_begin() const
{
  return dirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the dirichletbcs_ vector
//*******************************************************************|************************************************************//
std::vector< const dolfin::DirichletBC* >::iterator SystemBucket::dirichletbcs_end()
{
  return dirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the dirichletbcs_ vector
//*******************************************************************|************************************************************//
std::vector< const dolfin::DirichletBC* >::const_iterator SystemBucket::dirichletbcs_end() const
{
  return dirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the points_ vector
//*******************************************************************|************************************************************//
std::vector<ReferencePoints_ptr>::iterator SystemBucket::points_begin()
{
  return points_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the points_ vector
//*******************************************************************|************************************************************//
std::vector<ReferencePoints_ptr>::const_iterator SystemBucket::points_begin() const
{
  return points_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the points_ vector
//*******************************************************************|************************************************************//
std::vector<ReferencePoints_ptr>::iterator SystemBucket::points_end()
{
  return points_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the points_ vector
//*******************************************************************|************************************************************//
std::vector<ReferencePoints_ptr>::const_iterator SystemBucket::points_end() const
{
  return points_.end();
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
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the fields in the system
//*******************************************************************|************************************************************//
const std::string SystemBucket::fields_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionBucket_const_it f_it = fields_.begin();              // loop over the fields
                                    f_it != fields_.end(); f_it++ )
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
  for ( FunctionBucket_const_it f_it = coeffs_.begin();              // loop over the coefficients
                                  f_it != coeffs_.end(); f_it++ )
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
  for ( SolverBucket_const_it s_it = solvers_.begin();               // loop over the solvers
                                s_it != solvers_.end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// checkpoint the system
//*******************************************************************|************************************************************//
void SystemBucket::checkpoint()
{
  if (function())
  {

    std::stringstream buffer;                                        // optionpath buffer

    buffer.str(""); buffer << (*bucket()).output_basename() << "_" 
                           << name() << "_" 
                           << (*bucket()).checkpoint_count() << ".xml";
    dolfin::File sysfile(buffer.str());
    sysfile << *function();

  }

  for (FunctionBucket_it f_it = fields_begin();                      // if there's no function then there should be no fields
                         f_it != fields_end(); f_it++)               // so this is a bit redundant outside the above if statement
  {
    (*(*f_it).second).checkpoint();
  }

  checkpoint_options_();

}

//*******************************************************************|************************************************************//
// given a map from components to field initial condition expressions initialize the system initial condition expression
//*******************************************************************|************************************************************//
void SystemBucket::collect_ics_(const uint &components, const std::map< std::size_t, Expression_ptr > &icexpressions)
{
  const std::size_t nfields = icexpressions.size();
  if (nfields==1)                                                    // single field
  {
    const std::size_t rank = (*(*icexpressions.begin()).second).value_rank();
    if (rank==0)                                                     // scalar
    {
      icexpression_.reset(new InitialConditionExpression(icexpressions));
    }
    else if (rank==1)                                                // vector
    {
      icexpression_.reset(new InitialConditionExpression(components, icexpressions));
    }
    else if (rank==2)                                                // tensor
    {
      std::vector<std::size_t> value_shape(2, 0);
      for (uint i=0; i<rank; i++)
      {
        value_shape[i] = (*(*icexpressions.begin()).second).value_dimension(i);
      }
      icexpression_.reset(new InitialConditionExpression(value_shape, icexpressions));
    }
    else
    {
      dolfin::error("Unknown rank in collect_ics_");
    }
  }
  else                                                               // vectors are the general case for a mixed function
  {
    icexpression_.reset(new InitialConditionExpression(components, icexpressions));
  }
}

//*******************************************************************|************************************************************//
// initialize the system with a combined initial condition (calls eval)
//*******************************************************************|************************************************************//
void SystemBucket::apply_ic_()
{
  if (icexpression_)
  {
    (*oldfunction_).interpolate(*icexpression_);                     // interpolate the initial condition onto the old function
  }
  else if (icfile_)
  {
    (*icfile()) >> (*oldfunction_);
  }
  else
  {
    dolfin::error("Unknown way of applying initial condition.");
  }
  (*(*iteratedfunction_).vector()) = (*(*oldfunction_).vector());    // set the iterated function vector to the old function vector
  (*(*function_).vector()) = (*(*oldfunction_).vector());            // set the function vector to the old function vector
}

//*******************************************************************|************************************************************//
// apply the vector of system boundary conditions to the system function vectors to ensure consisten initial and boundary conditions
//*******************************************************************|************************************************************//
void SystemBucket::apply_dirichletbc_()
{
  for (int_FunctionBucket_const_it f_it = orderedfields_begin();     // loop over all the fields
                                f_it != orderedfields_end(); f_it++)
  {
    for (int_DirichletBC_const_it                                    // loop over all the bcs
         b_it = (*(*f_it).second).ordereddirichletbcs_begin(); 
         b_it != (*(*f_it).second).ordereddirichletbcs_end(); b_it++)
    {
      (*(*b_it).second).apply((*(*oldfunction_).vector()));
      (*(*b_it).second).apply((*(*iteratedfunction_).vector()));
      (*(*b_it).second).apply((*(*function_).vector()));
    }
  }
}

//*******************************************************************|************************************************************//
// apply the vector of system reference points to the system function vectors to ensure consistent initial and boundary conditions
//*******************************************************************|************************************************************//
void SystemBucket::apply_referencepoints_()
{
  for (int_FunctionBucket_const_it f_it = orderedfields_begin();     // loop over all the fields
                                f_it != orderedfields_end(); f_it++)
  {
    for (ReferencePoints_const_it                                    // loop over all the points
         p_it = (*(*f_it).second).points_begin(); 
         p_it != (*(*f_it).second).points_end(); p_it++)
    {
      (*(*p_it).second).apply((*(*oldfunction_).vector()));
      (*(*p_it).second).apply((*(*iteratedfunction_).vector()));
      (*(*p_it).second).apply((*(*function_).vector()));
    }
  }
}

//*******************************************************************|************************************************************//
// attach all coefficients, first to the functionals then to the solver forms
//*******************************************************************|************************************************************//
void SystemBucket::attach_all_coeffs_()
{
  attach_function_coeffs_(fields_begin(), fields_end());             // attach functions to field functionals
  attach_function_coeffs_(coeffs_begin(), coeffs_end());             // attach functions to coefficients functionals
  attach_solver_coeffs_(solvers_begin(), solvers_end());             // attach functions to the solver forms
}

//*******************************************************************|************************************************************//
// loop between the function bucket iterators attaching coefficients to the functionals of those function buckets
//*******************************************************************|************************************************************//
void SystemBucket::attach_function_coeffs_(FunctionBucket_it f_begin, 
                                             FunctionBucket_it f_end)
{
  for (FunctionBucket_it f_it = f_begin; f_it != f_end; f_it++)      // loop over the function buckets
  {
    (*(*f_it).second).attach_functional_coeffs();                    // attach coefficients to the functionals of this function
                                                                     // bucket
  }
}

//*******************************************************************|************************************************************//
// loop between the solver bucket iterators attaching coefficients to the forms of those solver buckets
//*******************************************************************|************************************************************//
void SystemBucket::attach_solver_coeffs_(SolverBucket_it s_begin, 
                                              SolverBucket_it s_end)
{

  for (SolverBucket_it s_it = s_begin; s_it != s_end; s_it++)        // loop over the solver buckets
  {
    (*(*s_it).second).attach_form_coeffs();                          // attach coefficients to the forms of this solver bucket
  }
}

//*******************************************************************|************************************************************//
// virtual checkpointing of options
//*******************************************************************|************************************************************//
void SystemBucket::checkpoint_options_()
{
  dolfin::error("Failed to find virtual function checkpoint_options_.");
}

//*******************************************************************|************************************************************//
// empty the data structures in the system bucket
//*******************************************************************|************************************************************//
void SystemBucket::empty_()
{
  fields_.clear();
  orderedfields_.clear();
  coeffs_.clear();
  orderedcoeffs_.clear();
  solvers_.clear();
  orderedsolvers_.clear();
}

