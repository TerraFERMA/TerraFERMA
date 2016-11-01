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
#include "FunctionBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include "BucketDolfinBase.h"
#include "DolfinPETScBase.h"
#include "BucketPETScBase.h"
#include "Logger.h"
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
  PetscErrorCode perr;                                               // petsc error code

  for (std::vector<IS>::iterator is = component_is_.begin(); 
                                 is != component_is_.end(); is++)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = ISDestroy(&(*is)); petsc_err(perr);                            // destroy the IS, necessary?
    #else
    perr = ISDestroy(*is); petsc_err(perr);                             // destroy the IS, necessary?
    #endif
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a generic function, which one depends on the time
//*******************************************************************|************************************************************//
const GenericFunction_ptr FunctionBucket::genericfunction_ptr(const double_ptr time) const
{
  if (time == (*(*system()).bucket()).current_time_ptr())
  {
    assert(iteratedfunction_);
    return iteratedfunction_;
  }
  else if (time == (*(*system()).bucket()).old_time_ptr())
  {
    assert(oldfunction_);
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
      assert(oldfunction_);
      return oldfunction_;
    }
  }
  else
  {
    tf_err("Unknown time pointer when returning function in genericfunction_ptr(time).", "FunctionBucket: %s, SystemBucket: %s", 
           name_.c_str(), (*system_).name().c_str());
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a generic function, which one depends on the function_type requested
//*******************************************************************|************************************************************//
const GenericFunction_ptr FunctionBucket::genericfunction_ptr(const std::string &function_type) const
{
  if (function_type=="iterated")
  {
    assert(iteratedfunction_);
    return iteratedfunction_;
  }
  else if (function_type=="old")
  {
    assert(oldfunction_);
    return oldfunction_;
  }
  else if (function_type=="residual")
  {
    assert(residualfunction_);
    return residualfunction_;
  }
  else if (function_type=="change")
  {
    assert(changefunction_);
    return changefunction_;
  }
  else if (function_type=="snesupdate")
  {
    assert(snesupdatefunction_);
    return snesupdatefunction_;
  }
  else
  {
    tf_err("Unknown function type when returning function in genericfunction_ptr(function_type).", 
           "function_type: %s, FunctionBucket name: %s, SystemBucket name: %s", 
           function_type.c_str(), name_.c_str(), (*system_).name().c_str());
  }
}

//*******************************************************************|************************************************************//
// cached the values vector (primarily for diagnostic purposes)
//*******************************************************************|************************************************************//
void FunctionBucket::cachevector(const std::string &function_type)
{
  if (!cachedvector_ || cachedvectortype_!=function_type)
  {
    clearcachedvector();
    cachedvector_ = basevector(function_type);
    cachedvectortype_ = function_type;
  }
}

//*******************************************************************|************************************************************//
// destroy the cached vector
//*******************************************************************|************************************************************//
void FunctionBucket::clearcachedvector()
{
  cachedvector_ = NULL;
  cachedvectortype_ = "no_cached_vector";
}

//*******************************************************************|************************************************************//
// return the base vector (could be the system vector for example) of values describing the function for the given function_type
//*******************************************************************|************************************************************//
const_PETScVector_ptr FunctionBucket::basevector(const std::string &function_type) const
{
  GenericFunction_ptr u = genericfunction_ptr(function_type);
  Mesh_ptr mesh = (*system_).mesh();
  const_PETScVector_ptr fv;

  const_Function_ptr uf = std::dynamic_pointer_cast<const dolfin::Function>(u);
  if (uf)
  {
    fv = std::dynamic_pointer_cast<const dolfin::PETScVector>((*uf).vector());
  }
  else
  {
    std::vector< double > values;
    (*u).compute_vertex_values(values, *mesh);

    std::size_t offset = 
                dolfin::MPI::global_offset((*mesh).mpi_comm(),
                                           values.size(),
                                           true);

    PETScVector_ptr tv( new dolfin::PETScVector() );
    (*tv).init((*mesh).mpi_comm(), std::make_pair(offset, offset+values.size()));
    (*tv).set_local(values);
    fv = std::const_pointer_cast<const dolfin::PETScVector>(tv);
  }
  return fv;
}

//*******************************************************************|************************************************************//
// return a vector of values describing the function for the given function_type
//*******************************************************************|************************************************************//
dolfin::PETScVector FunctionBucket::vector(const std::string &function_type, const int &component) const
{
  std::vector<int> components(1, component);
  return vector(function_type, &components);
}

//*******************************************************************|************************************************************//
// return a vector of values describing the function for the given function_type
//*******************************************************************|************************************************************//
dolfin::PETScVector FunctionBucket::vector(const std::string &function_type, const std::vector<int>* components) const
{
  PetscErrorCode perr;                                               // petsc error code
  Mesh_ptr mesh = (*system_).mesh();

  const_PETScVector_ptr fv;
  if (cachedvector_ && cachedvectortype_==function_type)
  {
    fv = cachedvector_;
  }
  else
  {
    fv = basevector(function_type);
  }

  IS is = components_is(components);
  PetscInt size;
  perr = ISGetLocalSize(is, &size);
  std::size_t offset = 
              dolfin::MPI::global_offset((*mesh).mpi_comm(),
                                         size, true);

  dolfin::PETScVector sv;
  sv.init((*mesh).mpi_comm(), std::make_pair(offset, offset+size));

  VecScatter scatter;
  perr = VecScatterCreate((*fv).vec(), is, 
                          sv.vec(), PETSC_NULL, 
                          &scatter);
  perr = VecScatterBegin(scatter, 
                         (*fv).vec(), sv.vec(), 
                         INSERT_VALUES, SCATTER_FORWARD);
  perr = VecScatterEnd(scatter,
                       (*fv).vec(), sv.vec(),
                       INSERT_VALUES, SCATTER_FORWARD);

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1            // necessary or taken care of when object leaves scope?
  perr = VecScatterDestroy(&scatter);
  perr = ISDestroy(&is);
  #else
  perr = VecScatterDestroy(scatter);
  perr = ISDestroy(is);
  #endif

  return sv;
}

//*******************************************************************|************************************************************//
// return an IS describing the requested component
// NOTE: IS should be destroyed after use
//*******************************************************************|************************************************************//
IS FunctionBucket::component_is(const int &component) const
{
  std::vector<int> components(1, component);
  return components_is(&components);
}

//*******************************************************************|************************************************************//
// return an IS describing the requested components
// NOTE: IS should be destroyed after use
//*******************************************************************|************************************************************//
IS FunctionBucket::components_is(const std::vector<int>* components) const
{
  Mesh_ptr mesh = (*system_).mesh();
  PetscErrorCode perr;                                               // petsc error code

  std::vector<int> subcomponents;
  if (components)
  {
    subcomponents = *components;
  }
  else
  {
    subcomponents.resize(size());
    std::iota(subcomponents.begin(), subcomponents.end(), 0);
  }

  IS islist[subcomponents.size()];
  for (uint i = 0; i < subcomponents.size(); i++)
  {
    islist[i] = component_is_[subcomponents[i]];
  }
  IS isall;
  perr = ISConcatenate((*mesh).mpi_comm(), subcomponents.size(), islist, &isall);

  return isall;
}

//*******************************************************************|************************************************************//
// return the maximum of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::max(const std::string &function_type, const uint component) const
{
  std::vector<int> components(1, component);
  return max(function_type, &components);
}

//*******************************************************************|************************************************************//
// return the maximum of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::max(const std::string &function_type, const std::vector<int>* components) const
{
  dolfin::PETScVector pvector = vector(function_type, components);
  return pvector.max();
}

//*******************************************************************|************************************************************//
// return the minimum of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::min(const std::string &function_type, const uint component) const
{
  std::vector<int> components(1, component);
  return min(function_type, &components);
}

//*******************************************************************|************************************************************//
// return the minimum of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::min(const std::string &function_type, const std::vector<int>* components) const
{
  dolfin::PETScVector pvector = vector(function_type, components);
  return pvector.min();
}

//*******************************************************************|************************************************************//
// return the norm of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::norm(const std::string &function_type, const std::string &norm_type, const uint component) const
{
  std::vector<int> components(1, component);
  return norm(function_type, norm_type, &components);
}

//*******************************************************************|************************************************************//
// return the norm of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::norm(const std::string &function_type, const std::string &norm_type, const std::vector<int>* components) const
{
  dolfin::PETScVector pvector = vector(function_type, components);
  return pvector.norm(norm_type);
}

//*******************************************************************|************************************************************//
// return the change in a component of the function bucket
//*******************************************************************|************************************************************//
double FunctionBucket::change(const uint component)
{
  std::vector<int> components(1, component);
  return change(&components);
}

//*******************************************************************|************************************************************//
// return the change in this function over a timestep (only valid for fields and only valid after system changefunction has been
// updated)
//*******************************************************************|************************************************************//
double FunctionBucket::change(const std::vector<int>* components)
{
  if(change_calculated_)                                             // just testing pointer is associated
  {
    assert(change_);
    if(!*change_calculated_)
    {
      assert(changefunction_);

      *change_ = norm("change", change_normtype_, components)/
                 norm("iterated", change_normtype_, components);

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
// return the dimension of the function
//*******************************************************************|************************************************************//
const std::size_t FunctionBucket::dimension(const std::size_t &i) const
{
  if (i >= shape_.size())
  {
    tf_err("Illegal dimension axis requested.", "Axis: %d, Rank: %d, FunctionBucket name: %s, SystemBucket name: %s",
           i, shape_.size(), name_.c_str(), (*system_).name().c_str());
  }
  return shape_[i];
}

//*******************************************************************|************************************************************//
// return if the function is symmetric or not (only true for fields)
//*******************************************************************|************************************************************//
const bool FunctionBucket::symmetric() const
{
  bool symmetric = false;
  if (functionspace_ && rank()==2)
  {
    symmetric = (*(*functionspace_).element()).num_sub_elements() != size();
  }
  return symmetric;
}

//*******************************************************************|************************************************************//
// refresh this functionbucket if it "needs" it
//*******************************************************************|************************************************************//
void FunctionBucket::refresh(const bool force)
{
  if (functiontype_==FUNCTIONBUCKET_FIELD)
  {
    std::vector<int> locations;
    locations.push_back(SOLVE_TIMELOOP);
    locations.push_back(SOLVE_DIAGNOSTICS);
    bool solved = (*system()).solve(locations, force);
  }
  else if (functiontype_==FUNCTIONBUCKET_COEFF)
  {
    if (coefficientfunction_)
    {
      (*std::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*coefficientfunction_);
    }
    if (constantfunctional_)
    {
      double value = dolfin::assemble(*constantfunctional_);
      *std::dynamic_pointer_cast< dolfin::Constant >(function_) = value;
    }
  }
  else
  {
    tf_err("Unknown FunctionBucket type.", "functiontype_ = %d", functiontype_);
  }
}

//*******************************************************************|************************************************************//
// update the potentially time dependent functions
//*******************************************************************|************************************************************//
void FunctionBucket::update()
{
  if (coefficientfunction_)
  {
    (*(*std::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).vector()) = 
      (*(*std::dynamic_pointer_cast< dolfin::Function >(function_)).vector());
                                                                     // update the oldfunction to the new function value
  }

  if (constantfunctional_)
  {
    *std::dynamic_pointer_cast< dolfin::Constant >(oldfunction_) = 
        double(*std::dynamic_pointer_cast< dolfin::Constant >(function_));
  }
}

//*******************************************************************|************************************************************//
// update the calculated flags
//*******************************************************************|************************************************************//
void FunctionBucket::resetcalculated()
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
    (*std::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*coefficientfunction_);
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
    *std::dynamic_pointer_cast< dolfin::Constant >(function_) = value;
  }
}

//*******************************************************************|************************************************************//
// cap the values in the vector associated with this functionspace
//*******************************************************************|************************************************************//
void FunctionBucket::postprocess_values()
{
  assert(functiontype_ == FUNCTIONBUCKET_FIELD);

  const std::size_t lsize = size();

  bool cap = false;
  bool zero = false;
  for(uint i = 0; i < lsize; i++)
  {
    if (upper_cap_[i] || lower_cap_[i])
    {
      cap = true;
    }
    if (zeropoints_[i])
    {
      zero = true;
    }
  }

  if (!(cap||zero))
  {
    return;
  }

  PetscErrorCode perr;                                               // petsc error code
  PETScVector_ptr sysvec = std::dynamic_pointer_cast<dolfin::PETScVector>((*(*system()).iteratedfunction()).vector());

  for (uint i = 0; i < lsize; i++)
  {
    if (!(zeropoints_[i]||upper_cap_[i]||lower_cap_[i]))
    {
      continue;
    }

    double value = 0.0;
    if (zeropoints_[i])
    {
      std::vector< Array_double_ptr > values;
      (*zeropoints_[i]).eval(values, *function(), (*system_).mesh());
      double lvalue = 0.0;
      if (values.size()>0)
      {
        assert(values.size()==1);
        lvalue = (*values[0])[i];
      }
      std::size_t nvals = dolfin::MPI::sum((*(*system_).mesh()).mpi_comm(),  // we calculate this just in case more than 1 process
                                   values.size());                   // found the same value
      assert(nvals>0);
      value = dolfin::MPI::sum((*(*system_).mesh()).mpi_comm(),
                               lvalue)/nvals;
    }

    PetscInt np;
    const PetscInt *pindices;
    perr = ISGetLocalSize(component_is_[i], &np);
    petsc_err(perr);
    perr = ISGetIndices(component_is_[i], &pindices);
    petsc_err(perr);

    for (uint j = 0; j < np; j++)
    {
      double current_value;
      (*sysvec).get(&current_value, 1, &pindices[j]);

      if (upper_cap_[i])
      {
        if (current_value > *upper_cap_[i])
        {
          (*sysvec).setitem(pindices[j], *upper_cap_[i]);
        }
      }

      if (lower_cap_[i])
      {
        if (current_value < *lower_cap_[i])
        {
          (*sysvec).setitem(pindices[j], *lower_cap_[i]);
        }
      }

      if (zeropoints_[i])
      {
        (*sysvec).setitem(pindices[j], current_value-value);
      }
    }

    perr = ISRestoreIndices(component_is_[i], &pindices);
    petsc_err(perr);

  }

  (*sysvec).apply("insert");
  *(*(*system_).function()).vector() = *sysvec;
}

//*******************************************************************|************************************************************//
// loop over the functionals in this function bucket and attach the coefficients they request using the parent bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::attach_form_coeffs()
{
  if (system_)
  {
    if (constantfunctional_)
    {
      (*(*system_).bucket()).attach_coeffs(constantfunctional_);
    }
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc expression in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_bcfunction(GenericFunction_ptr bcfunction, const std::string &name)
{
  GenericFunction_it g_it = bcfunctions_.find(name);                      // check if the name already exists
  if (g_it != bcfunctions_.end())
  {
    tf_err("BCFunction already exists in function.", "BCExpression name: %s, Function name: %s, System name: %s", 
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    bcfunctions_[name] = bcfunction;                                 // if not, then insert the expression pointer into the maps
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a bc expression form from the function bucket data maps
//*******************************************************************|************************************************************//
GenericFunction_ptr FunctionBucket::fetch_bcfunction(const std::string &name)
{
  GenericFunction_it g_it = bcfunctions_.find(name);                    // check if the name already exists
  if (g_it == bcfunctions_.end())
  {
    tf_err("BCFunction does not exist in function.", "BCExpression name: %s, Function name: %s, System name: %s", 
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*g_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a bc in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_dirichletbc(DirichletBC_ptr bc, const std::string &name)
{
  DirichletBC_hash_it bc_it = dirichletbcs_.get<om_key_hash>().find(name);                                  // check if the name already exists
  if (bc_it != dirichletbcs_.get<om_key_hash>().end())
  {
    tf_err("DirichletBC already exists in function.", "DirichletBC name: %s, Function name: %s, System name: %s", 
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    dirichletbcs_.insert(om_item<const std::string, DirichletBC_ptr>(name, bc));                                        // if not, register the bc
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::dirichletbcs_begin()
{
  return dirichletbcs_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::dirichletbcs_begin() const
{
  return dirichletbcs_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_it FunctionBucket::dirichletbcs_end()
{
  return dirichletbcs_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the dirichletbcs_ map
//*******************************************************************|************************************************************//
DirichletBC_const_it FunctionBucket::dirichletbcs_end() const
{
  return dirichletbcs_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a reference point in the function bucket data maps
//*******************************************************************|************************************************************//
void FunctionBucket::register_referencepoint(ReferencePoint_ptr point, const std::string &name)
{
  ReferencePoint_hash_it p_it = referencepoints_.get<om_key_hash>().find(name);             // check if the name already exists
  if (p_it != referencepoints_.get<om_key_hash>().end())
  {
    tf_err("ReferencePoint already exists in function.", "ReferencePoint name: %s, Function name: %s, System name: %s", 
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    referencepoints_.insert(om_item<const std::string,ReferencePoint_ptr>(name, point));             // if not, register the bc
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the referencepoints_ map
//*******************************************************************|************************************************************//
ReferencePoint_it FunctionBucket::referencepoints_begin()
{
  return referencepoints_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the referencepoints_ map
//*******************************************************************|************************************************************//
ReferencePoint_const_it FunctionBucket::referencepoints_begin() const
{
  return referencepoints_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the referencepoints_ map
//*******************************************************************|************************************************************//
ReferencePoint_it FunctionBucket::referencepoints_end()
{
  return referencepoints_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the referencepoints_ map
//*******************************************************************|************************************************************//
ReferencePoint_const_it FunctionBucket::referencepoints_end() const
{
  return referencepoints_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// include this function in visualization output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_visualization() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_visualization.");
  return false;
}

//*******************************************************************|************************************************************//
// include the residual of this function in visualization output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_residual_in_visualization() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_residual_in_visualization.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in diagnostic output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_statistics() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_statistics.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in steadystate output and checking
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_steadystate() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_steadystate.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in detectors output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionBucket::include_in_detectors() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_detectors.");
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
  return s.str();
}

//*******************************************************************|************************************************************//
// fill the component is's
//*******************************************************************|************************************************************//
void FunctionBucket::fill_is_()
{
  const std::size_t lsize = size();
  Mesh_ptr mesh = (*system_).mesh();

  component_is_.resize(lsize);

  if (functionspace())                                               // field and function coefficients should end up here
  {
    const std::size_t lrank = rank();
    if (lrank==0||lrank==1)
    {
      for (int i = 0; i < lsize; i++)
      {
        std::vector<std::size_t> indices = local_functionspace_dofs(functionspace(), i);
        component_is_[i] = convert_vector_to_is((*mesh).mpi_comm(), 
                                                indices);
      }
    }
    else if (lrank==2)
    {
      const std::size_t dim0 = dimension(0);
      const std::size_t dim1 = dimension(1);
      const bool lsymmetric = symmetric();
      for (int i = 0; i < dim0; i++)
      {
        for (int j = 0; j < dim1; j++)
        {
          std::size_t k = i*dim1 + j;
          std::size_t kf = k;
          if (lsymmetric)
          {
            if (j >= i)
            {
              kf = i*dim1 + j - (i*(i+1))/2;
            }
            else
            {
              kf = j*dim1 + i - (j*(j+1))/2;
            }
          }

          std::vector<std::size_t> indices = local_functionspace_dofs(functionspace(), kf);
          component_is_[k] = convert_vector_to_is((*mesh).mpi_comm(), 
                                                  indices);
        }
      }
    }
    else
    {
      tf_err("Unknown rank in fill_is_.", "Rank: %d", lrank);
    }
  }
  else                                                               // only expression coefficients and constants
  {                                                                  // should end up here, in which case we're just
                                                                     // indexing into a vector of vertex values
      const std::size_t nv = (*mesh).num_vertices();
      std::size_t offset = 
                  dolfin::MPI::global_offset((*mesh).mpi_comm(),
                                             lsize*nv,
                                             true);
      for (int i = 0; i < lsize; i++)
      {
        std::vector<std::size_t> indices(nv);
        std::iota(indices.begin(), indices.end(), offset + i*(*mesh).num_vertices());
        component_is_[i] = convert_vector_to_is((*mesh).mpi_comm(), 
                                                indices);
      }
    
  }
  
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
  tf_err("Failed to find virtual function checkpoint_options_.", "Need to implement a checkpointing method.");
}

