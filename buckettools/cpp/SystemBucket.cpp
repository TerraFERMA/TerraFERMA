
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
  }
}

//*******************************************************************|************************************************************//
// update the timelevels of the system function
// this is only done for the system function at the moment as this is the only time dependent function
// field functions are subfunctions and therefore require deep copies to be made anyway while coefficients are currently not
// time dependent
//*******************************************************************|************************************************************//
void SystemBucket::update()
{
  (*oldfunction_).vector() = (*function_).vector();
}

//*******************************************************************|************************************************************//
// attach coefficients to forms and functionals then initialize matrices described by this system's forms
//*******************************************************************|************************************************************//
void SystemBucket::attach_and_initialize()
{
  dolfin::info("Attaching coeffs for system %s", name().c_str());
  attach_all_coeffs_();                                              // attach the coefficients to form and functionals

  for (SolverBucket_it s_it = solvers_begin(); s_it != solvers_end();// loop over the solver buckets
                                                              s_it++)
  {
    (*(*s_it).second).initialize_matrices();                         // perform a preassembly of all the matrices to set up
                                                                     // sparsities etc.
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_field(FunctionBucket_ptr field, std::string name)
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
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a field function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_field(std::string name)
{
  FunctionBucket_it f_it = fields_.find(name);                       // check if name already exists
  if (f_it == fields_.end())
  {
    dolfin::error("Field named \"%s\" does not exists in system.",   // if it doesn't, issue an error
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
const FunctionBucket_ptr SystemBucket::fetch_field(std::string name) const
{
  FunctionBucket_const_it f_it = fields_.find(name);                 // check if name already exists
  if (f_it == fields_.end())
  {
    dolfin::error("Field named \"%s\" does not exists in system.",   // if it doesn't, issue an error
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
// register a (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_coeff(FunctionBucket_ptr coeff, std::string name)
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
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a coefficient function bucket in the system bucket data maps
//*******************************************************************|************************************************************//
FunctionBucket_ptr SystemBucket::fetch_coeff(std::string name)
{
  FunctionBucket_it f_it = coeffs_.find(name);                       // check if the name already exists
  if (f_it == coeffs_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Coefficient named \"%s\" does not exists in system.", 
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
// register a (boost shared) pointer to a solver bucket in the system bucket data maps
//*******************************************************************|************************************************************//
void SystemBucket::register_solver(SolverBucket_ptr solver, std::string name)
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
void SystemBucket::output()
{
  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end(); 
                                                              f_it++)
  {
    (*(*f_it).second).output();
  }
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
const std::string SystemBucket::fields_str(int indent) const
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
const std::string SystemBucket::coeffs_str(int indent) const
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
const std::string SystemBucket::solvers_str(int indent) const
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
  coeffs_.clear();
  solvers_.clear();
  orderedsolvers_.clear();
}

