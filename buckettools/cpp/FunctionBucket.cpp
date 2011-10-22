
#include "BoostTypes.h"
#include "FunctionBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
FunctionBucket::FunctionBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
FunctionBucket::FunctionBucket(SystemBucket* system) : system_(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
FunctionBucket::~FunctionBucket()
{
  empty_();                                                          // empty the function data maps
}

//*******************************************************************|************************************************************//
// return the change in this function over a timestep (only valid for fields and only valid after system changefunction has been
// updated)
//*******************************************************************|************************************************************//
const double FunctionBucket::change()
{
  if(change_calculated_)
  {
    assert(change_);
    if(!*change_calculated_)
    {
      assert(changefunction_);

      dolfin::Function changefunc =                                    // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >(changefunction());
      *change_ = changefunc.vector().norm(change_normtype_);

      *change_calculated_=true;
    }
    return *change_;
  }
  else
  {
    return 0.0;
  }
}

//*******************************************************************|************************************************************//
// reset the change boolean
//*******************************************************************|************************************************************//
void FunctionBucket::resetchange()
{
  if(change_calculated_)
  {
    *change_calculated_=false;
  }
}

//*******************************************************************|************************************************************//
// update the potentially time dependent functions
//*******************************************************************|************************************************************//
void FunctionBucket::update_timedependent()
{
  if (coefficientfunction_)
  {
    (*boost::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).vector() = 
      (*boost::dynamic_pointer_cast< dolfin::Function >(function_)).vector();
                                                                     // update the oldfunction to the new function value
    (*boost::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*coefficientfunction_);
  }
}

//*******************************************************************|************************************************************//
// update the potentially nonlinear functions
//*******************************************************************|************************************************************//
void FunctionBucket::update_nonlinear()
{
  if (constantfunctional_)
  {
    double value = dolfin::assemble(*constantfunctional_);
    *boost::dynamic_pointer_cast< dolfin::Constant >(function_) = value;
  }
}

//*******************************************************************|************************************************************//
// loop over the functionals in this function bucket and attach the coefficients they request using the parent bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::attach_functional_coeffs()
{
  if (system_)
  {
    if (constantfunctional_)
    {
      (*(*system_).bucket()).attach_coeffs(constantfunctional_);
    }

    (*(*system_).bucket()).attach_coeffs(functionals_begin(),
                                                  functionals_end());
  }
}

//*******************************************************************|************************************************************//
// make a partial copy of the provided function bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void FunctionBucket::copy_diagnostics(FunctionBucket_ptr &function, SystemBucket_ptr &system) const
{

  if(!function)
  {
    function.reset( new FunctionBucket(&(*system)) );
  }

  (*function).name_ = name_;
  (*function).index_ = index_;

  (*function).functionspace_ = functionspace_;

  (*function).function_ = function_;
  (*function).iteratedfunction_ = iteratedfunction_;
  (*function).oldfunction_ = oldfunction_;

  (*function).changefunction_ = changefunction_;
  (*function).change_ = change_;
  (*function).change_calculated_ = change_calculated_;
  (*function).change_normtype_ = change_normtype_;

  (*function).functionals_ = functionals_;
  (*function).bcexpressions_ = bcexpressions_;
  (*function).bcs_ = bcs_;

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a functional form in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_functional(Form_ptr functional, const std::string &name)
{
  Form_it f_it = functionals_.find(name);                            // check if the name already exists
  if (f_it != functionals_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
            "Functional named \"%s\" already exists in function.", 
                                                      name.c_str());
  }
  else
  {
    functionals_[name] = functional;                                 // if not, insert it into the functionals_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a functional form from the function bucket data maps
//*******************************************************************|************************************************************//
Form_ptr FunctionBucket::fetch_functional(const std::string &name)
{
  Form_it f_it = functionals_.find(name);                            // check if the name already exists
  if (f_it == functionals_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Functional named \"%s\" does not exist in function.", 
                                                      name.c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the functionals_ map
//*******************************************************************|************************************************************//
Form_it FunctionBucket::functionals_begin()
{
  return functionals_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the functionals_ map
//*******************************************************************|************************************************************//
Form_const_it FunctionBucket::functionals_begin() const
{
  return functionals_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the functionals_ map
//*******************************************************************|************************************************************//
Form_it FunctionBucket::functionals_end()
{
  return functionals_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the functionals_ map
//*******************************************************************|************************************************************//
Form_const_it FunctionBucket::functionals_end() const
{
  return functionals_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc expression in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_bcexpression(Expression_ptr bcexpression, const std::string &name)
{
  Expression_it e_it = bcexpressions_.find(name);                    // check if the name already exists
  if (e_it != bcexpressions_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
            "BCExpression named \"%s\" already exists in function.", 
                                                      name.c_str());
  }
  else
  {
    bcexpressions_[name] = bcexpression;                             // if not, then insert the expression pointer into the maps
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_bc(BoundaryCondition_ptr bc, const std::string &name)
{
  BoundaryCondition_it bc_it = bcs_.find(name);                      // check if the name already exists
  if (bc_it != bcs_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
        "BoundaryCondition named \"%s\" already exists in function.", 
                                                    name.c_str());
  }
  else
  {
    bcs_[name] = bc;                                                 // if not, register the bc
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the bcs_ map
//*******************************************************************|************************************************************//
BoundaryCondition_it FunctionBucket::bcs_begin()
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the bcs_ map
//*******************************************************************|************************************************************//
BoundaryCondition_const_it FunctionBucket::bcs_begin() const
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the bcs_ map
//*******************************************************************|************************************************************//
BoundaryCondition_it FunctionBucket::bcs_end()
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the bcs_ map
//*******************************************************************|************************************************************//
BoundaryCondition_const_it FunctionBucket::bcs_end() const
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// output the current contents of the function to a pvd file (if associated)
//*******************************************************************|************************************************************//
void FunctionBucket::output(const bool &write_vis)
{
  if (pvdfile_ && write_vis)                                         // check a pvd file is associated
  {
    *pvdfile_ << *boost::dynamic_pointer_cast< dolfin::Function >(function());
  }
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the function bucket
//*******************************************************************|************************************************************//
const std::string FunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the functionals in the function bucket
//*******************************************************************|************************************************************//
const std::string FunctionBucket::functionals_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = functionals_.begin(); f_it != functionals_.end(); f_it++ )
  {
    s << indentation << "Functional " << (*f_it).first  << std::endl;
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// include this function in visualization output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_visualization() const
{
  dolfin::error("Failed to find virtual function include_in_visualization.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in diagnostic output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_statistics() const
{
  dolfin::error("Failed to find virtual function include_in_statistics.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in steadystate output and checking
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_steadystate() const
{
  dolfin::error("Failed to find virtual function include_in_steadystate.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in detectors output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_detectors() const
{
  dolfin::error("Failed to find virtual function include_in_steadystate.");
  return false;
}

//*******************************************************************|************************************************************//
// empty the data structures in the function bucket
//*******************************************************************|************************************************************//
void FunctionBucket::empty_()
{
  functionals_.clear();
  bcexpressions_.clear();
  bcs_.clear();
}

