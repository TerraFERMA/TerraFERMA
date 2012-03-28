
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
    boost::unordered_set<uint> dofs = (*(*functionspace()).dofmap()).dofs();

    boost::unordered_set<uint>::const_iterator dof;
    for (dof = dofs.begin(); dof != dofs.end(); dof++)
    {
      if(upper_cap_)
      {
        if ((*(*(*system()).function()).vector())[*dof] > *upper_cap_)
        {
          (*(*(*system()).function()).vector()).setitem(*dof, *upper_cap_);
        }
      }

      if(lower_cap_)
      {
        if ((*(*(*system()).function()).vector())[*dof] < *lower_cap_)
        {
          (*(*(*system()).function()).vector()).setitem(*dof, *lower_cap_);
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

  (*function).functionals_ = functionals_;
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
    orderedbcs_[(int) bcs_.size()] = bc;                             // also insert it in the order it was registered in the 
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
// return an iterator to the beginning of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_BoundaryCondition_it FunctionBucket::orderedbcs_begin()
{
  return orderedbcs_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_BoundaryCondition_const_it FunctionBucket::orderedbcs_begin() const
{
  return orderedbcs_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_BoundaryCondition_it FunctionBucket::orderedbcs_end()
{
  return orderedbcs_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedbcs_ map
//*******************************************************************|************************************************************//
int_BoundaryCondition_const_it FunctionBucket::orderedbcs_end() const
{
  return orderedbcs_.end();
}

//*******************************************************************|************************************************************//
// output the current contents of the function to a pvd file (if associated)
//*******************************************************************|************************************************************//
void FunctionBucket::output(const bool &write_vis)
{
  if (write_vis)
  {
    if (pvdfile_)                                                    // check a pvd file is associated
    {
      *pvdfile_ << std::make_pair(&(*boost::dynamic_pointer_cast< dolfin::Function >(function())),
                                  (*(*system()).bucket()).current_time());
    }
    if (respvdfile_)                                                 // check a residual pvd file is associated
    {
      *respvdfile_ << std::make_pair(&(*boost::dynamic_pointer_cast< dolfin::Function >(residualfunction())),
                                     (*(*system()).bucket()).current_time());
    }
  }
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
// include this function in detectors output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_detectors() const
{
  dolfin::error("Failed to find virtual function include_in_steadystate.");
  return false;
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

