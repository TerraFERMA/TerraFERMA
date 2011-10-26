
#include "BoostTypes.h"
#include "InitialConditionExpression.h"
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

  (*oldfunction_).vector() = (*function_).vector();                  // update the oldfunction to the new function value

                                                                     // fields share a vector with the system function so no need to
                                                                     // update them...

  update_nonlinear();                                                // update potentially nonlinear coefficients

  resetchange();                                                     // reset the change booleans in the system and fields

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
  }
  return maxchange;
}

//*******************************************************************|************************************************************//
// update the change function
//*******************************************************************|************************************************************//
void SystemBucket::updatechange()
{
  if (!*change_calculated_)
  {
    assert(changefunction_);

    (*changefunction_).vector() = (*function_).vector();             // before updating the oldfunction to the new values
    (*changefunction_).vector() -= (*oldfunction_).vector();         // update the change in the fields over this timesteps

    *change_calculated_ = true;
  }
}

//*******************************************************************|************************************************************//
// reset the change function
//*******************************************************************|************************************************************//
void SystemBucket::resetchange()
{
  *change_calculated_ = false;
  for (FunctionBucket_it f_it = fields_begin(); 
                                  f_it != fields_end(); f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())
    {
      (*(*f_it).second).resetchange();
    }
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
// return an iterator to the beginning of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector<BoundaryCondition_ptr>::iterator SystemBucket::bcs_begin()
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector<BoundaryCondition_ptr>::const_iterator SystemBucket::bcs_begin() const
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector<BoundaryCondition_ptr>::iterator SystemBucket::bcs_end()
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the bcs_ vector
//*******************************************************************|************************************************************//
std::vector<BoundaryCondition_ptr>::const_iterator SystemBucket::bcs_end() const
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// loop over the fields outputting pvd diagnostics for all the fields in this system
//*******************************************************************|************************************************************//
void SystemBucket::output(const bool &write_vis)
{
  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    (*(*f_it).second).output(write_vis);
  }
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
    include = (*(*f_it).second).include_in_visualization();
    if (include)
    {
      break;
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
// loop over all the fields in this system, collecting their bcs into a single vector of system bcs
//*******************************************************************|************************************************************//
void SystemBucket::collect_bcs_()
{
  for (FunctionBucket_const_it f_it = fields_begin();                // loop over all the fields
                                        f_it != fields_end(); f_it++)
  {
    for (BoundaryCondition_const_it                                  // loop over all the bcs
          b_it = (*(*f_it).second).bcs_begin(); 
          b_it != (*(*f_it).second).bcs_end(); b_it++)
    {
      bcs_.push_back((*b_it).second);                                // add the bcs to a std vector
    }
  }
}

//*******************************************************************|************************************************************//
// given a map from components to field initial condition expressions initialize the system initial condition expression
//*******************************************************************|************************************************************//
void SystemBucket::collect_ics_(const uint &component, 
              const std::map< uint, Expression_ptr > &icexpressions)
{
  if (component==1)
  {
    icexpression_.reset(new InitialConditionExpression(icexpressions));// the system function is scalar so set up a scalar ic expression
  }
  else
  {                                                                  // multiple components so set up a multi-component ic
                                                                     // expression
    icexpression_.reset( new InitialConditionExpression(component, icexpressions));
  }
}

//*******************************************************************|************************************************************//
// initialize the system with a combined initial condition (calls eval)
//*******************************************************************|************************************************************//
void SystemBucket::apply_ic_()
{
  (*oldfunction_).interpolate(*icexpression_);                       // interpolate the initial condition onto the old function
  (*iteratedfunction_).vector() = (*oldfunction_).vector();          // set the iterated function vector to the old function vector
  (*function_).vector() = (*oldfunction_).vector();                  // set the function vector to the old function vector
}

//*******************************************************************|************************************************************//
// apply the vector of system boundary conditions to the system function vectors to ensure consisten initial and boundary conditions
//*******************************************************************|************************************************************//
void SystemBucket::apply_bc_()
{
  for (std::vector<BoundaryCondition_ptr>::const_iterator            // loop over the vector of bcs
                            bc = bcs_begin(); bc != bcs_end(); bc++)
  {
    (*(*bc)).apply((*oldfunction_).vector());
    (*(*bc)).apply((*iteratedfunction_).vector());
    (*(*bc)).apply((*function_).vector());
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

