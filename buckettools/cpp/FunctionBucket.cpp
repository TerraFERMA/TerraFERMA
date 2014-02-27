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


#include "SystemBucket.h"
#include "BoostTypes.h"
#include "FunctionBucket.h"
#include "Bucket.h"
#include "BucketDolfinBase.h"
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
// return a pointer to a generic function, which one depends on the time
//*******************************************************************|************************************************************//
const GenericFunction_ptr FunctionBucket::genericfunction_ptr(const double_ptr time) const
{
  if (time == (*(*system()).bucket()).current_time_ptr())
  {
    return iteratedfunction_;
  }
  else if (time == (*(*system()).bucket()).old_time_ptr())
  {
    return oldfunction_;
  }
  else if (time == (*(*system()).bucket()).start_time_ptr())
  {
    if (icexpression_)
    {
      return icexpression_;                                          // if possible, return the initial condition expression
    }
    else                                                             // but coefficients don't have icexpressions so let's assume
    {                                                                // we're only evaluating this as an initial condition...
      assert((*time-DOLFIN_EPS)<=(*(*system()).bucket()).old_time());// check we're at the old time and this is valid (FIXME: not a
                                                                     // sufficient check as this should only be called in the 
                                                                     // initialization routine, where this is guaranteed anyway!)
      return oldfunction_;
    }
  }
  else
  {
    dolfin::error("Unknown time pointer when returning function in genericfunction(time).");
  }
}

//*******************************************************************|************************************************************//
// return the maximum of the function bucket at the given time
//*******************************************************************|************************************************************//
const double FunctionBucket::max(const double_ptr time, const int *index0, const int *index1) const
{
  const GenericFunction_ptr function = genericfunction_ptr(time);
  double maxvalue;

  if (  functiontype_==FUNCTIONBUCKET_FIELD || 
      ((functiontype_==FUNCTIONBUCKET_COEFF)&&(type()=="Function")) )
  {
    dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
      *boost::dynamic_pointer_cast< const dolfin::Function >(function);

    if (index0 && index1)
    {
      assert(func.value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      dolfin::Function funccomp = func[*index0][*index1];
      maxvalue = (*funccomp.vector()).max();
    }
    else if (index0)
    {
      assert(func.value_rank()==1);
      assert(*index0 < (*function).value_size());

      dolfin::Function funccomp = func[*index0];
      maxvalue = (*funccomp.vector()).max();
    }
    else
    {
      maxvalue = (*func.vector()).max();
    }
  }
  else
  {
    Mesh_ptr mesh = (*system()).mesh();
    const int nvert = (*mesh).num_vertices();
    std::vector< double > values;
    (*function).compute_vertex_values(values, *mesh);

    if (index0 && index1)
    {
      assert((*function).value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      maxvalue = *std::max_element(&values[(2*(*index1)+(*index0))*nvert], 
                  &values[(2*(*index1)+(*index0)+1)*nvert]);  // maximum for requested component
    }
    else if (index0)
    {
      assert((*function).value_rank()==1);
      assert(*index0 < (*function).value_size());

      maxvalue = *std::max_element(&values[(*index0)*nvert], 
                  &values[((*index0)+1)*nvert]);              // maximum for requested component
    }
    else
    {
      maxvalue = *std::max_element(&values[0], &values[values.size()]);
    }
  }

  return maxvalue;
  
}

//*******************************************************************|************************************************************//
// return the minimum of the function bucket at the given time
//*******************************************************************|************************************************************//
const double FunctionBucket::min(const double_ptr time, const int *index0, const int *index1) const
{
  const GenericFunction_ptr function = genericfunction_ptr(time);
  double minvalue;

  if (  functiontype_==FUNCTIONBUCKET_FIELD || 
      ((functiontype_==FUNCTIONBUCKET_COEFF)&&(type()=="Function")) )
  {
    dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
      *boost::dynamic_pointer_cast< const dolfin::Function >(function);

    if (index0 && index1)
    {
      assert(func.value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      dolfin::Function funccomp = func[*index0][*index1];
      minvalue = (*funccomp.vector()).min();
    }
    else if (index0)
    {
      assert(func.value_rank()==1);
      assert(*index0 < (*function).value_size());

      dolfin::Function funccomp = func[*index0];
      minvalue = (*funccomp.vector()).min();
    }
    else
    {
      minvalue = (*func.vector()).min();
    }
  }
  else
  {
    Mesh_ptr mesh = (*system()).mesh();
    const int nvert = (*mesh).num_vertices();
    std::vector< double > values;
    (*function).compute_vertex_values(values, *mesh);

    if (index0 && index1)
    {
      assert((*function).value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      minvalue = *std::min_element(&values[(2*(*index1)+(*index0))*nvert], 
                  &values[(2*(*index1)+(*index0)+1)*nvert]);  // minimum for requested component
    }
    else if (index0)
    {
      assert((*function).value_rank()==1);
      assert(*index0 < (*function).value_size());

      minvalue = *std::min_element(&values[(*index0)*nvert], 
                  &values[((*index0)+1)*nvert]);              // minimum for requested component
    }
    else
    {
      minvalue = *std::min_element(&values[0], &values[values.size()]);
    }
  }

  return minvalue;
  
}

//*******************************************************************|************************************************************//
// return the infnorm of the function bucket at the given time
//*******************************************************************|************************************************************//
const double FunctionBucket::infnorm(const double_ptr time, const int *index0, const int *index1) const
{
  const GenericFunction_ptr function = genericfunction_ptr(time);
  double normvalue;

  if (  functiontype_==FUNCTIONBUCKET_FIELD || 
      ((functiontype_==FUNCTIONBUCKET_COEFF)&&(type()=="Function")) )
  {
    dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
      *boost::dynamic_pointer_cast< const dolfin::Function >(function);

    if (index0 && index1)
    {
      assert(func.value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      dolfin::Function funccomp = func[*index0][*index1];
      normvalue = (*funccomp.vector()).norm("linf");
    }
    else if (index0)
    {
      assert(func.value_rank()==1);
      assert(*index0 < (*function).value_size());

      dolfin::Function funccomp = func[*index0];
      normvalue = (*funccomp.vector()).norm("linf");
    }
    else
    {
      normvalue = (*func.vector()).norm("linf");
    }
  }
  else
  {
    Mesh_ptr mesh = (*system()).mesh();
    const int nvert = (*mesh).num_vertices();
    std::vector< double > values;
    (*function).compute_vertex_values(values, *mesh);

    if (index0 && index1)
    {
      assert((*function).value_rank()==2);
      assert(*index0 < (*function).value_dimension(0));
      assert(*index1 < (*function).value_dimension(1));

      normvalue = std::abs(*std::max_element(&values[(2*(*index1)+(*index0))*nvert], 
                           &values[(2*(*index1)+(*index0)+1)*nvert], abslessthan));
    }
    else if (index0)
    {
      assert((*function).value_rank()==1);
      assert(*index0 < (*function).value_size());

      normvalue = std::abs(*std::max_element(&values[(*index0)*nvert], 
                           &values[((*index0)+1)*nvert], abslessthan));
    }
    else
    {
      normvalue = std::abs(*std::max_element(&values[0], 
                           &values[nvert], abslessthan));
    }
  }

  return normvalue;
  
}

//*******************************************************************|************************************************************//
// return the maximum of the function at the current time
//*******************************************************************|************************************************************//
const double FunctionBucket::functionmax() const
{
  return max((*(*system()).bucket()).current_time_ptr());
}

//*******************************************************************|************************************************************//
// return the minimum of the function at the current time
//*******************************************************************|************************************************************//
const double FunctionBucket::functionmin() const
{
  return min((*(*system()).bucket()).current_time_ptr());
}

//*******************************************************************|************************************************************//
// return the infnorm of the function at the current time
//*******************************************************************|************************************************************//
const double FunctionBucket::functioninfnorm() const
{
  return infnorm((*(*system()).bucket()).current_time_ptr());
}

//*******************************************************************|************************************************************//
// return the size of the function
//*******************************************************************|************************************************************//
const std::size_t FunctionBucket::size() const
{
  std::size_t size = 1;
  for (std::vector<std::size_t>::const_iterator i=shape_.begin(); i != shape_.end(); i++)
  {
    size*=(*i);
  }
  return size;
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
      *change_ = (*changefunc.vector()).norm(change_normtype_);

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
// return the change in the value of the given functional over a timestep
//*******************************************************************|************************************************************//
const double FunctionBucket::functionalchange(Form_const_it f_it)
{
  return std::abs(functionalvalue(f_it)-oldfunctionalvalue(f_it));
}

//*******************************************************************|************************************************************//
// return the value of the given functional (calculating it if necessary)
//*******************************************************************|************************************************************//
const double FunctionBucket::functionalvalue(Form_const_it f_it)
{
  bool_ptr calculated = functional_calculated_[(*f_it).first];
  if(*calculated)
  {
    return *functional_values_[(*f_it).first];
  }
  else
  {
    double value = dolfin::assemble((*(*f_it).second));              // assemble the functional
    *functional_values_[(*f_it).first] = value;
    *calculated = true;

    return value;
  }
}

//*******************************************************************|************************************************************//
// return the old stored value of the given functional
//*******************************************************************|************************************************************//
const double FunctionBucket::oldfunctionalvalue(Form_const_it f_it)
{
  return *oldfunctional_values_[(*f_it).first];
}

//*******************************************************************|************************************************************//
// reset the calculated booleans
//*******************************************************************|************************************************************//
void FunctionBucket::resetcalculated()
{
  if(change_calculated_)
  {
    *change_calculated_=false;
  }

  for (bool_ptr_it b_it = functional_calculated_.begin(); 
                       b_it != functional_calculated_.end(); b_it++)
  {
    *(*b_it).second = false;
  }
}

//*******************************************************************|************************************************************//
// refresh this functionbucket if it "needs" it
//*******************************************************************|************************************************************//
void FunctionBucket::refresh(const bool force)
{
  if (functiontype_==FUNCTIONBUCKET_FIELD)
  {
    if( (!(*system()).solved()) || force )                           // solve the field's system if it hasn't already
    {                                                                // been solved for this timestep or we're forcing it
      (*system()).solve();
    }
  }
  else if (functiontype_==FUNCTIONBUCKET_COEFF)
  {
    if (coefficientfunction_)
    {
      (*boost::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*coefficientfunction_);
    }
    if (constantfunctional_)
    {
      double value = dolfin::assemble(*constantfunctional_);
      *boost::dynamic_pointer_cast< dolfin::Constant >(function_) = value;
    }
  }
  else
  {
    dolfin::error("Unknown FunctionBucket type.");
  }
}

//*******************************************************************|************************************************************//
// update the potentially time dependent functions
//*******************************************************************|************************************************************//
void FunctionBucket::update()
{
  if (coefficientfunction_)
  {
    (*(*boost::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).vector()) = 
      (*(*boost::dynamic_pointer_cast< dolfin::Function >(function_)).vector());
                                                                     // update the oldfunction to the new function value
  }

  if (constantfunctional_)
  {
    *boost::dynamic_pointer_cast< dolfin::Constant >(oldfunction_) = 
        double(*boost::dynamic_pointer_cast< dolfin::Constant >(function_));
  }

  for (Form_const_it f_it = functionals_begin(); f_it != functionals_end(); f_it++)
  {
    *oldfunctional_values_[(*f_it).first] = *functional_values_[(*f_it).first];
  }
}

//*******************************************************************|************************************************************//
// update the potentially time dependent functions
//*******************************************************************|************************************************************//
void FunctionBucket::update_timedependent()
{
  if (coefficientfunction_)
  {
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
// cap the values in the vector associated with this functionspace
//*******************************************************************|************************************************************//
void FunctionBucket::cap_values()
{
  if (lower_cap_ || upper_cap_)
  {

    std::pair<uint, uint> ownership_range =                          // the parallel ownership range of the system functionspace
                      (*(*functionspace()).dofmap()).ownership_range();
    for (uint i = ownership_range.first; i < ownership_range.second; i++)
    {
      if(upper_cap_)
      {
        if ((*(*(*system()).function()).vector())[i] > *upper_cap_)
        {
          (*(*(*system()).function()).vector()).setitem(i, *upper_cap_);
        }
      }

      if(lower_cap_)
      {
        if ((*(*(*system()).function()).vector())[i] < *lower_cap_)
        {
          (*(*(*system()).function()).vector()).setitem(i, *lower_cap_);
        }
      }

      (*(*(*system()).function()).vector()).apply("insert");
    }
 
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

  (*function).residualfunction_ = residualfunction_;
  (*function).snesupdatefunction_ = snesupdatefunction_;

  (*function).functionals_ = functionals_;
  (*function).functional_values_ = functional_values_;
  (*function).oldfunctional_values_ = oldfunctional_values_;
  (*function).functional_calculated_ = functional_calculated_;

  (*function).bcexpressions_ = bcexpressions_;
  (*function).bcs_ = bcs_;
  (*function).points_ = points_;

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
    double_ptr value;
    value.reset( new double(0.0) );
    functional_values_[name]      = value;
    value.reset( new double(0.0) );
    oldfunctional_values_[name]   = value;
    bool_ptr calculated;
    calculated.reset( new bool(false) );
    functional_calculated_[name]  = calculated;
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
// return a (boost shared) pointer to a bc expression form from the function bucket data maps
//*******************************************************************|************************************************************//
Expression_ptr FunctionBucket::fetch_bcexpression(const std::string &name)
{
  Expression_it e_it = bcexpressions_.find(name);                    // check if the name already exists
  if (e_it == bcexpressions_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "BCExpression named \"%s\" does not exist in function.", 
                                                      name.c_str());
  }
  else
  {
    return (*e_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_bc(DirichletBC_ptr bc, const std::string &name)
{
  DirichletBC_it bc_it = bcs_.find(name);                      // check if the name already exists
  if (bc_it != bcs_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
        "BoundaryCondition named \"%s\" already exists in function.", 
                                                    name.c_str());
  }
  else
  {
    bcs_[name] = bc;                                                 // if not, register the bc
    orderedbcs_[(int) bcs_.size()] = bc;                             // also insert it in the order it was registered in the 
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the bcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::bcs_begin()
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the bcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::bcs_begin() const
{
  return bcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the bcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::bcs_end()
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the bcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::bcs_end() const
{
  return bcs_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_it FunctionBucket::orderedbcs_begin()
{
  return orderedbcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_const_it FunctionBucket::orderedbcs_begin() const
{
  return orderedbcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_it FunctionBucket::orderedbcs_end()
{
  return orderedbcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_const_it FunctionBucket::orderedbcs_end() const
{
  return orderedbcs_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_dirichletbc(DirichletBC_ptr bc, const std::string &name)
{
  DirichletBC_it bc_it;
  
  bc_it = dirichletbcs_.find(name);                                  // check if the name already exists
  if (bc_it != dirichletbcs_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
        "Dirichlet BoundaryCondition named \"%s\" already exists in function.", 
                                                    name.c_str());
  }
  else
  {
    dirichletbcs_[name] = bc;                                        // if not, register the bc
    ordereddirichletbcs_[(int) dirichletbcs_.size()] = bc;           // also insert it in the order it was registered in the 
  }

  register_bc(bc, name);                                             // register the dirichlet bc in the normal maps too
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::dirichletbcs_begin()
{
  return dirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::dirichletbcs_begin() const
{
  return dirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::dirichletbcs_end()
{
  return dirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::dirichletbcs_end() const
{
  return dirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the ordereddirichletbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_it FunctionBucket::ordereddirichletbcs_begin()
{
  return ordereddirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the ordereddirichletbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_const_it FunctionBucket::ordereddirichletbcs_begin() const
{
  return ordereddirichletbcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the ordereddirichletbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_it FunctionBucket::ordereddirichletbcs_end()
{
  return ordereddirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the ordereddirichletbcs_ map
//*******************************************************************|************************************************************//
int_DirichletBC_const_it FunctionBucket::ordereddirichletbcs_end() const
{
  return ordereddirichletbcs_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a reference point in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_point(ReferencePoints_ptr point, const std::string &name)
{
  ReferencePoints_it p_it = points_.find(name);                       // check if the name already exists
  if (p_it != points_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
        "ReferencePoints named \"%s\" already exists in function.", 
                                                    name.c_str());
  }
  else
  {
    points_[name] = point;                                           // if not, register the bc
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the points_ map
//*******************************************************************|************************************************************//
ReferencePoints_it FunctionBucket::points_begin()
{
  return points_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the points_ map
//*******************************************************************|************************************************************//
ReferencePoints_const_it FunctionBucket::points_begin() const
{
  return points_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the points_ map
//*******************************************************************|************************************************************//
ReferencePoints_it FunctionBucket::points_end()
{
  return points_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the points_ map
//*******************************************************************|************************************************************//
ReferencePoints_const_it FunctionBucket::points_end() const
{
  return points_.end();
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
// include the residual of this function in visualization output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_residual_in_visualization() const
{
  dolfin::error("Failed to find virtual function include_residual_in_visualization.");
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
// include this a functional of this function with the given name in steadystate output and checking
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_functional_in_steadystate(const std::string &name) const
{
  dolfin::error("Failed to find virtual function include_functional_in_steadystate.");
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
// checkpoint the functionbucket
//*******************************************************************|************************************************************//
void FunctionBucket::checkpoint()
{
  checkpoint_options_();
}

//*******************************************************************|************************************************************//
// virtual checkpointing of options
//*******************************************************************|************************************************************//
void FunctionBucket::checkpoint_options_()
{
  dolfin::error("Failed to find virtual function checkpoint_options_.");
}

//*******************************************************************|************************************************************//
// empty the data structures in the function bucket
//*******************************************************************|************************************************************//
void FunctionBucket::empty_()
{
  functionals_.clear();
  bcexpressions_.clear();
  bcs_.clear();
  orderedbcs_.clear();
}

