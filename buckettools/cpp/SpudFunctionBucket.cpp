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


#include "PythonExpression.h"
#include "BoostTypes.h"
#include "SpudFunctionBucket.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemFunctionalsWrapper.h"
#include "SystemExpressionsWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "RegionsExpression.h"
#include "SemiLagrangianExpression.h"
#include "PointDetectors.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudFunctionBucket::SpudFunctionBucket(const std::string &optionpath, 
                                            SystemBucket* system) : 
                                            optionpath_(optionpath), 
                                            FunctionBucket(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudFunctionBucket::~SpudFunctionBucket()
{
}

//*******************************************************************|************************************************************//
// fill the function bucket data structures assuming the buckettools schema and that this is a field
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_field(const uint &index)
{
  functiontype_ = FUNCTIONBUCKET_FIELD;                              // this is a field

  fill_base_(index);                                                 // fill the base data about the function bucket
                                                                     // (could be called from constructor as it's common regardless
                                                                     // of functionbucket type)

                                                                     // check we think this is a field (necessary but not
                                                                     // sufficient)
  assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field::"+name()));

  allocate_field_();                                                 // allocate as a field

  fill_is_();                                                        // fill information about the component index sets
}

//*******************************************************************|************************************************************//
// fill the function bucket data structures assuming the buckettools schema and that this is a coefficient
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_coeff(const uint &index)
{
  functiontype_ = FUNCTIONBUCKET_COEFF;                              // this is a coefficient

  fill_base_(index);                                                 // fill the base data about the function bucket
                                                                     // (called from constructor instead?)

  if(system_)                                                        // can't do this if this coefficient lives outside a system
  { 
                                                                     // check we think this is a coefficient (not a foolproof test)
    assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/coefficient::"+name()));
  }

  allocate_coeff_expression_();                                      // allocate as a coefficient (doesn't do much for coefficient
                                                                     // functions)

}

//*******************************************************************|************************************************************//
// allocate a coefficient that is of type function - left until almost last as we need a functionspace for it
//*******************************************************************|************************************************************//
void SpudFunctionBucket::allocate_coeff_function()
{

  std::stringstream buffer;

  if (type_=="Function")                                             // only do this for coefficient functions
  {
    if(!(*(*system_).bucket()).contains_coefficientspace(uflsymbol()))// check we have a coefficient space for this function's ufl
    {                                                                // symbol
      tf_err("Unable to allocate coefficient as no matching functionspace found.",
             "Coefficient %s declared as type::Function but its ufl_symbol was not found in any forms or functionals.",
             name().c_str());
    }

    functionspace_ =                                                 // grab the functionspace for this coefficient from the bucket
          (*(*system_).bucket()).fetch_coefficientspace(uflsymbol());// data maps
    outputfunctionspace_ = functionspace_;

    function_.reset( new dolfin::Function(functionspace_) );        // allocate the function on this functionspace
    oldfunction_.reset( new dolfin::Function(functionspace_) );     // allocate the old function on this functionspace
    iteratedfunction_ = function_;                                   // just point this at the function

                                                                     // can't initialize this yet (it may depend on other
                                                                     // expressions through a user defined cpp expression) so just
                                                                     // zero if for now and we'll intialize it later, once we're
                                                                     // allowed to call eval for the first time
    (*(*std::dynamic_pointer_cast< dolfin::Function >(function_)).vector()).zero();
    (*(*std::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).vector()).zero();

    buffer.str(""); buffer << (*system_).name() << "::" << name();     // rename the function as SystemName::CoefficientName
    (*function_).rename(buffer.str(), buffer.str());

    buffer.str(""); buffer << (*system_).name() << "::Old" << name();  // rename the old function as SystemName::OldCoefficientName
    (*oldfunction_).rename(buffer.str(), buffer.str());

    fill_is_();                                                       // we finally have a functionspace so can set up is's

  }

}

//*******************************************************************|************************************************************//
// allocate bcs
//*******************************************************************|************************************************************//
void SpudFunctionBucket::allocate_bcs()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath()                             // find out how many strong bcs there are (weak bcs handled in
                                << "/type/rank/boundary_condition";  // forms)
  int nbcs = Spud::option_count(buffer.str());
  if (nbcs > 0)                                                      // if we have any...
  {
    for (uint i = 0; i < nbcs; i++)                                  // loop over the bcs
    {
      buffer.str(""); buffer << optionpath()
                    << "/type/rank/boundary_condition[" << i << "]";
      fill_bc_(buffer.str());                                        // and fill in details about the bc
    }
  }

                                                                     // NOTE: this isn't strictly necessary here but we do it
                                                                     // for consistency as these are basically bcs too
  buffer.str(""); buffer << optionpath()                             // find out how many reference points bcs there are
                                << "/type/rank/reference_point";
  int nrpoints = Spud::option_count(buffer.str());
  if (nrpoints > 0)                                                  // if we have any...
  {
    for (uint i = 0; i < nrpoints; i++)                              // loop over the points
    {
      buffer.str(""); buffer << optionpath()
                    << "/type/rank/reference_point[" << i << "]";
      fill_reference_point_(buffer.str());                           // and fill in details about the reference point
    }
  }
    
}

//*******************************************************************|************************************************************//
// initialize the function bucket data structures (mostly expressions) assuming the buckettools schema and that this is a field
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_field()
{
  std::stringstream buffer;                                          // optionpath buffer

  buffer.str(""); buffer << optionpath()                             // find out how many strong bcs there are (weak bcs handled in
                                << "/type/rank/boundary_condition";  // forms)
  int nbcs = Spud::option_count(buffer.str());
  for (uint i = 0; i < nbcs; i++)                                    // loop over the bcs
  {
    buffer.str(""); buffer << optionpath()
                  << "/type/rank/boundary_condition[" << i << "]";
    initialize_bc_(buffer.str());                                    // and initialize the expressions describing the bc
  }
    
  buffer.str(""); buffer << optionpath()
                   << "/type/rank/initial_condition::WholeMesh/file";
  if(!Spud::have_option(buffer.str()))
  {
    buffer.str(""); buffer << optionpath()                           // find out how many initial conditions we have
                                    << "/type/rank/initial_condition";
    if(Spud::have_option(buffer.str()))
    {
      initialize_expression_over_regions_(icexpression_, buffer.str());// intialize the expression
    }
  }
  
}

//*******************************************************************|************************************************************//
// initialize the function bucket data structures (mostly expressions) assuming the buckettools schema and that this is a coefficient
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_coeff_expression()
{
  std::stringstream buffer;                                          // optionpath buffer

  if (type_=="Expression" || type_=="Constant")                      // only proceed for expressions and constants
  {

    buffer.str(""); buffer << optionpath() << "/type/rank/value";
    initialize_expression_over_regions_(
        std::dynamic_pointer_cast< dolfin::Expression >(function_), 
                                                   buffer.str());    // iteratedfunction_ points at this too so doesn't need
                                                                     // initializing
    initialize_expression_over_regions_(
        std::dynamic_pointer_cast< dolfin::Expression >(oldfunction_), 
                                                   buffer.str());
    
  }
}

//*******************************************************************|************************************************************//
// initialize the function bucket data structures (mostly expressions) assuming the buckettools schema and that this is a coefficient
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_coeff_function()
{
  std::stringstream buffer;                                          // optionpath buffer

  if (type_=="Function")                                             // only proceed for coefficient functions
  {
    buffer.str(""); buffer << optionpath() << "/type/rank/value";

    bool time_dependent;
    Expression_ptr tmpexpression = allocate_expression_over_regions_(
                                   buffer.str(), 
                                   (*(*system_).bucket()).current_time_ptr(), 
                                   &time_dependent);

    initialize_expression_over_regions_(tmpexpression, buffer.str());

    (*std::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*tmpexpression);
    (*std::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).interpolate(*tmpexpression);
                                                                     // iteratedfunction_ points at function_
    if (time_dependent)
    {
      coefficientfunction_ = tmpexpression;                          // we'll need this again
    }

  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in visualization output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_visualization() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_visualization");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the residual of this function bucket should be included in visualization output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_residual_in_visualization() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_residual_in_visualization");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in diagnostic output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_statistics() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_statistics");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in steady state check and output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_steadystate() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_steady_state");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in detectors output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_detectors() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_detectors");
}

//*******************************************************************|************************************************************//
// return a string describing the contents of this function bucket
//*******************************************************************|************************************************************//
const std::string SpudFunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << " (" 
                                << optionpath() << ")" << std::endl;
  return s.str();
}

//*******************************************************************|************************************************************//
// fill the function bucket base data assuming the buckettools schema (common for fields and coefficients)
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_base_(const uint &index)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  index_ = index;                                                    // the field or coefficient index (most relevant for fields)

  buffer.str(""); buffer << optionpath() << "/name";                 // field or coefficient name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/ufl_symbol";           // ufl symbol
  serr = Spud::get_option(buffer.str(), uflsymbol_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/name";            // field or coefficient type (function, expression, constant)
  serr = Spud::get_option(buffer.str(), type_);                      // (as string)
  spud_err(buffer.str(), serr);

  std::string strrank;
  buffer.str(""); buffer << optionpath() << "/type/rank/name";       // field or coefficient rank (as string)
  serr = Spud::get_option(buffer.str(), strrank); 
  spud_err(buffer.str(), serr);

  shape_.resize(0);                                                  // start by assuming scalar

  if (strrank=="Vector")                                             // if it's a vector we find the shape
  {                                                                  // 2 ways depending on if it's a constant
    int size_int;                                                    // or an expression
    if (type_=="Constant")
    {
      buffer.str(""); buffer << optionpath()
                                      << "/type/rank/value/constant";
      std::vector< int > shape_int;
      serr = Spud::get_option_shape(buffer.str(), shape_int);
      spud_err(buffer.str(), serr);
      size_int = shape_int[0];
    }
    else
    {
      buffer.str(""); buffer << optionpath()                         // the size of the element (if a vector)
                                        << "/type/rank/element/size";// NOTE: defaults to the bucket geometry dimension
      serr = Spud::get_option(buffer.str(), size_int, (*(*system_).bucket()).dimension()); 
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << optionpath()
                                      << "/type/rank/value";
      int nvalues = Spud::option_count(buffer.str());
      for (uint i = 0; i < nvalues; i++)
      {
        buffer.str(""); buffer << optionpath()
                                        << "/type/rank/value[" << i << "]/constant";
        if (Spud::have_option(buffer.str()))
        {
          std::vector< int > constant_shape;
          serr = Spud::get_option_shape(buffer.str(), constant_shape);
          spud_err(buffer.str(), serr);
          assert(size_int==constant_shape[0]);                       // just checking!
        }
      }
    }
    shape_.resize(1);
    shape_[0] = size_int;
  }

  if (strrank=="Tensor")
  {
    std::vector< int > shape_int;
    if (type_=="Constant")
    {
      buffer.str(""); buffer << optionpath()
                                      << "/type/rank/value/constant";
      serr = Spud::get_option_shape(buffer.str(), shape_int); 
      spud_err(buffer.str(), serr);
    }
    else
    {
      buffer.str(""); buffer << optionpath()                         // the shape of the element (if a tensor)
                                       << "/type/rank/element/shape";// NOTE: defaults to the bucket geometry dimension square
      std::vector< int > default_shape(2, (*(*system_).bucket()).dimension());
      serr = Spud::get_option(buffer.str(), shape_int, default_shape); 
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << optionpath()
                                      << "/type/rank/value";
      int nvalues = Spud::option_count(buffer.str());
      for (uint i = 0; i < nvalues; i++)
      {
        buffer.str(""); buffer << optionpath()
                                        << "/type/rank/value[" << i << "]/constant";
        if (Spud::have_option(buffer.str()))
        {
          std::vector< int > constant_shape;
          serr = Spud::get_option_shape(buffer.str(), constant_shape);
          spud_err(buffer.str(), serr);
          assert((shape_int[0]==constant_shape[0])&&(shape_int[1]==constant_shape[1]));// just checking!
        }
      }
    }
    shape_.resize(2);
    shape_[0] = shape_int[0];
    shape_[1] = shape_int[1];
  }
 
}

//*******************************************************************|************************************************************//
// allocate a field assuming the base data has already been filled and using a buckettools spud schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::allocate_field_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  assert(type_=="Function");                                         // fields may only be functions

  int nfields =                                                      // how may fields in the system (i.e. mixed functionspace?)
      Spud::option_count((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field");
  if (nfields == 1)                                                  // no subfunctions (this field is identical to the system)
  {                                                                  // just grab references to:
    functionspace_ = (*system_).functionspace();                     // the functionspace
    outputfunctionspace_ = functionspace_;


    function_ = (*system_).function();                               // the function
    oldfunction_ = (*system_).oldfunction();                         // the old function
    iteratedfunction_ = (*system_).iteratedfunction();               // the iterated function

    if (include_in_steadystate())
    {
      changefunction_ = (*system_).changefunction();                 // the change in the function between timesteps
    }

    residualfunction_ = (*system_).residualfunction();

    if ((*system_).snesupdatefunction())
    {
      snesupdatefunction_ = (*system_).snesupdatefunction();
    }

  }
  else                                                               // yes, multiple fields in this system so we need to make
  {                                                                  // a subspace and subfunctions (dangerous, be careful!)
    functionspace_ = (*(*system_).functionspace())[index_];          // create a subspace
    outputfunctionspace_ = (*functionspace_).collapse();             // the collapsed functionspace for this function

    function_.reset( &(*(*system_).function())[index_],              // create a subfunction but don't delete it ever as it shares
                                              dolfin::NoDeleter() ); // information with the system function
    oldfunction_.reset( &(*(*system_).oldfunction())[index_],        // same for the old function
                                              dolfin::NoDeleter() );
    iteratedfunction_.reset( &(*(*system_).iteratedfunction())[index_], 
                                              dolfin::NoDeleter() ); // and the iterated function

    if (include_in_steadystate())
    {
      changefunction_.reset( &(*(*system_).changefunction())[index_],
                                              dolfin::NoDeleter() ); // and the change in the function between timesteps
    }

    residualfunction_.reset( &(*(*system_).residualfunction())[index_],
                                              dolfin::NoDeleter() );
    if ((*system_).snesupdatefunction())
    {
      snesupdatefunction_.reset( &(*(*system_).snesupdatefunction())[index_], 
                                              dolfin::NoDeleter() );
    }

  }

  tosystem_.reset( new dolfin::FunctionAssigner(functionspace_, 
                                                outputfunctionspace_) );
  fromsystem_.reset( new dolfin::FunctionAssigner(outputfunctionspace_, 
                                                  functionspace_) );

  buffer.str(""); buffer << (*system_).name() << "::" << name();     // rename the function as SystemName::FieldName
  (*function_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << (*system_).name() << "::Old" << name();  // rename the old function as SystemName::OldFieldName
  (*oldfunction_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << (*system_).name() << "::Iterated"        // rename the iterated function as SystemName::IteratedFieldName
                                                         << name();  
  (*iteratedfunction_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << (*system_).name() << "::Residual"        // rename the residual function as SystemName::ResidualFieldName
                                                         << name();  
  (*residualfunction_).rename(buffer.str(), buffer.str());

  if (include_in_steadystate())
  {
    buffer.str(""); buffer << (*system_).name() << "::TimestepChange"// rename the change function as SystemName::TimestepChangeFieldName
                                                           << name();  
    (*changefunction_).rename(buffer.str(), buffer.str());

    change_.reset( new double(0.0) );

    change_calculated_.reset( new bool(false) );

    buffer.str(""); buffer << optionpath()
                      << "/diagnostics/include_in_steady_state/norm";
    serr = Spud::get_option(buffer.str(), change_normtype_);
    spud_err(buffer.str(), serr);
  }

  if ((*system_).snesupdatefunction())
  {
    buffer.str(""); buffer << (*system_).name() << "::SNESUpdate"    // rename the snes update function as SystemName::SNESUpdateFieldName
                                                          << name();  
    (*snesupdatefunction_).rename(buffer.str(), buffer.str());
  }

  zeropoints_.resize(size());
  buffer.str(""); buffer << optionpath()                             // find out how many zero points bcs there are
                                << "/type/rank/zero_point";
  int nzpoints = Spud::option_count(buffer.str());
  if (nzpoints > 0)                                                  // if we have any...
  {
    for (uint i = 0; i < nzpoints; i++)                              // loop over the points
    {
      buffer.str(""); buffer << optionpath()
                    << "/type/rank/zero_point[" << i << "]";
      fill_zero_point_(buffer.str());                                // and fill in details about the zero point
    }
  }
    
  buffer.str(""); buffer << optionpath()
                   << "/type/rank/initial_condition::WholeMesh/file";
  if(Spud::have_option(buffer.str()))
  {
    serr = Spud::get_option(buffer.str(), icfilename_);
    spud_err(buffer.str(), serr);
  }
  else
  {
    buffer.str(""); buffer << optionpath()                           // set up initial condition optionpath
                                   << "/type/rank/initial_condition";// if this doesn't exist allocate_expression will
    icexpression_ = allocate_expression_over_regions_(buffer.str(),  // default to a 0 constant...
                          (*(*system_).bucket()).start_time_ptr());  // allocate the expression
    icfilename_ = "";
  }

  lower_cap_.resize(size());
  upper_cap_.resize(size());
  buffer.str(""); buffer << optionpath()                             // find out how many value caps there are
                                << "/type/rank/value_cap";
  int ncaps = Spud::option_count(buffer.str());
  if (ncaps > 0)                                                     // if we have any...
  {
    for (uint i = 0; i < ncaps; i++)                                 // loop over the caps
    {
      buffer.str(""); buffer << optionpath()
                    << "/type/rank/value_cap[" << i << "]";
      fill_cap_(buffer.str());                                       // and fill in details about the cap
    }
  }
    
}

//*******************************************************************|************************************************************//
// fill the details of a bc assuming the buckettools schema and given a set of edge ids in a mesh function
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_bc_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string bcname;                                                // bc name
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), bcname); 
  spud_err(buffer.str(), serr);
  
  std::vector<int> bcids;                                            // the boundary ids this bc applies to
  buffer.str(""); buffer << optionpath << "/boundary_ids";
  serr = Spud::get_option(buffer.str(), bcids); 
  spud_err(buffer.str(), serr);
  
  buffer.str(""); buffer << optionpath << "/sub_components";         // how many subcomponents exist?
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsubcomp; i++)                                // loop over the subcomponents
  {
    buffer.str(""); buffer << optionpath << "/sub_components[" 
                                                        << i << "]";
    fill_bc_component_(buffer.str(), bcname, bcids);                 // initialize this component
  }
}

//*******************************************************************|************************************************************//
// fill the details of a component of a bc on some boundary ids  assuming the buckettools schema and given a set of edge ids in a 
// mesh function
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_bc_component_(const std::string &optionpath,
                                            const std::string &bcname,
                                            const std::vector<int> &bcids)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;                                      // a buffer to assemble this bc's component name

  std::string compname;
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), compname);                   // get the name of these subcomponents
  spud_err(buffer.str(), serr);

  std::string bctype;
  buffer.str(""); buffer << optionpath << "/type/name";              // get the bc type
  serr = Spud::get_option(buffer.str(), bctype);
  spud_err(buffer.str(), serr);

  if(bctype=="Dirichlet")
  {
    namebuffer.str(""); namebuffer << bcname << compname;              // assemble a name for the bc expression
    buffer.str(""); buffer << optionpath << "/type::Dirichlet";      // initialize an expression describing the function
    GenericFunction_ptr bcfunc = allocate_expression_(buffer.str(), 
                     namebuffer.str(), 
                     (*(*system_).bucket()).current_time_ptr());     // on this boundary condition, assuming it is of type
    if (!bcfunc)                                                     // if nothing was returned then it's probably a 
    {                                                                // reference, which we handle separately
      bcfunc = take_function_reference_(buffer.str(), 
                    namebuffer.str(), 
                    (*(*system_).bucket()).current_time_ptr());
    }
    assert(bcfunc);

    register_bcfunction(bcfunc, namebuffer.str());                  // register the expression in the function bucket maps
       
    buffer.str(""); buffer << optionpath << "/components";           // do we have any components (i.e. could be scalar)?
    if (Spud::have_option(buffer.str()))
    {                                                                
      std::vector<int> subcompids;
      serr = Spud::get_option(buffer.str(), subcompids);             // get a list of the subcomponents 
      spud_err(buffer.str(), serr);
      
      for (std::vector<int>::const_iterator subcompid =              // loop over those subcomponents
                                            subcompids.begin(); 
                        subcompid < subcompids.end(); subcompid++)
      {

         FunctionSpace_ptr subfunctionspace =                        // grab the subspace from the field functionspace
                                      (*functionspace())[*subcompid];// happily dolfin caches this for us
         
         for (std::vector<int>::const_iterator bcid = bcids.begin(); // loop over the boundary ids
                                            bcid < bcids.end(); bcid++)
         {                                                           // create a new bc on each boundary id for this subcomponent
           DirichletBC_ptr bc(new dolfin::DirichletBC(subfunctionspace, bcfunc, (*system_).facetdomains(), *bcid));
           namebuffer.str(""); namebuffer << bcname << "::" 
                                      << *subcompid << "::" << *bcid;// assemble a name incorporating the boundary id
           register_dirichletbc(bc, namebuffer.str());               // register the bc in the function bucket
         }
         
      }
    }
    else                                                             // no components (scalar or using all components)
    {
      for (std::vector<int>::const_iterator bcid = bcids.begin();    // loop over the boundary ids
                                          bcid < bcids.end(); bcid++)
      {                                                              // create a bc on each boundary id for all components
        DirichletBC_ptr bc(new dolfin::DirichletBC(functionspace(), bcfunc, (*system_).facetdomains(), *bcid));
        namebuffer.str(""); namebuffer << bcname << "::" << *bcid;   // assemble a name for this bc incorporating the boundary id
        register_dirichletbc(bc, namebuffer.str());                  // register the bc in the function bucket
      }
    }
  }
  else
  {
    tf_err("Unknown bc type.", "bctype = ", bctype.c_str());
  }

}

//*******************************************************************|************************************************************//
// fill the details of a reference point assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_reference_point_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;                                      // a buffer to assemble this point's component name
  
  std::string pointname;                                             // point name
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), pointname); 
  spud_err(buffer.str(), serr);
  
  std::vector<double> coord;                                         // the coordinate of the reference point
  buffer.str(""); buffer << optionpath << "/coordinates";
  serr = Spud::get_option(buffer.str(), coord);
  spud_err(buffer.str(), serr);
  
  Constant_ptr value( new dolfin::Constant(0.0) );                   // works because we should always be scalar below

  const uint num_sub_elements = (*(*functionspace()).element()).num_sub_elements();
  if (num_sub_elements>0)
  {
    std::vector<int> subcompids;
    buffer.str(""); buffer << optionpath 
                           << "/sub_components/components";          // do subcomponents exist?
    if (Spud::have_option(buffer.str()))
    {
      serr = Spud::get_option(buffer.str(), subcompids);             // get a list of the subcomponents 
      spud_err(buffer.str(), serr);
    }
    else
    {
      subcompids.resize(size());                                     // we do this to ensure we're applying
      std::iota(subcompids.begin(), subcompids.end(), 0);            // ref point to scalar functionspace
    }
      
    for (std::vector<int>::const_iterator subcompid =                // loop over those subcomponents
                                          subcompids.begin(); 
                      subcompid < subcompids.end(); subcompid++)
    {

      FunctionSpace_ptr subfunctionspace =                           // grab the subspace from the field functionspace
                                    (*functionspace())[*subcompid];  // happily dolfin caches this for us
       
      ReferencePoint_ptr point(new ReferencePoint(coord, subfunctionspace, value));
      namebuffer.str(""); namebuffer << pointname << "::" 
                                 << *subcompid;                      // assemble a name incorporating the subcompid
      register_referencepoint(point, namebuffer.str());              // register the point in the function bucket
       
    }
  }
  else                                                               // no components (scalar or using all components)
  {                                                                  // this case is included because you can't index a scalar fs
    ReferencePoint_ptr point(new ReferencePoint(coord, functionspace(), value));
    register_referencepoint(point, pointname);                       // register the point in the function bucket
  }
}

//*******************************************************************|************************************************************//
// fill the details of a zero point assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_zero_point_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string pointname;                                             // point name
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), pointname); 
  spud_err(buffer.str(), serr);
  
  std::vector<double> coord;                                         // the coordinate of the reference point
  buffer.str(""); buffer << optionpath << "/coordinates";
  serr = Spud::get_option(buffer.str(), coord);
  spud_err(buffer.str(), serr);

  PointDetectors_ptr zeropoint( new PointDetectors(coord, pointname) );
  
  std::vector<int> subcompids;
  buffer.str(""); buffer << optionpath 
                         << "/sub_components/components";            // how many subcomponents exist?
  if (Spud::have_option(buffer.str()))
  {
    serr = Spud::get_option(buffer.str(), subcompids);               // get a list of the subcomponents 
    spud_err(buffer.str(), serr);
  }
  else
  {
    subcompids.resize(size());
    std::iota(subcompids.begin(), subcompids.end(), 0);
  }
    
  for (std::vector<int>::const_iterator subcompid =                // loop over those subcomponents
                                        subcompids.begin(); 
                    subcompid < subcompids.end(); subcompid++)
  {
    if (zeropoints_[*subcompid])
    {
      log(INFO, "WARNING: Multiple zero points applied to %s::%s (component %d).  Taking latest.",
                                (*system_).name().c_str(), name().c_str(), *subcompid);
    }
    zeropoints_[*subcompid] = zeropoint;
  }
}

//*******************************************************************|************************************************************//
// fill the details of a value cap assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_cap_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  double_ptr upper_cap;
  buffer.str(""); buffer << optionpath << "/upper_cap";
  if(Spud::have_option(buffer.str()))
  {
    upper_cap.reset( new double );
    serr = Spud::get_option(buffer.str(), *upper_cap);
    spud_err(buffer.str(), serr);
  }
  
  double_ptr lower_cap;
  buffer.str(""); buffer << optionpath << "/lower_cap";
  if(Spud::have_option(buffer.str()))
  {
    lower_cap.reset( new double );
    serr = Spud::get_option(buffer.str(), *lower_cap);
    spud_err(buffer.str(), serr);
  }
  
  std::vector<int> subcompids;
  buffer.str(""); buffer << optionpath 
                         << "/sub_components/components";            // how many subcomponents exist?
  if (Spud::have_option(buffer.str()))
  {
    serr = Spud::get_option(buffer.str(), subcompids);               // get a list of the subcomponents 
    spud_err(buffer.str(), serr);
  }
  else
  {
    subcompids.resize(size());
    std::iota(subcompids.begin(), subcompids.end(), 0);
  }
    
  for (std::vector<int>::const_iterator subcompid =                // loop over those subcomponents
                                        subcompids.begin(); 
                    subcompid < subcompids.end(); subcompid++)
  {
    assert(*subcompid<size());

    if (upper_cap)
    {
      if (upper_cap_[*subcompid])
      {
        log(INFO, "WARNING: Multiple upper caps applied to %s::%s (component %d).  Taking minimum.",
                                  (*system_).name().c_str(), name().c_str(), *subcompid);
        *upper_cap_[*subcompid] = std::min(*upper_cap_[*subcompid], *upper_cap);
      }
      else
      {
        upper_cap_[*subcompid] = upper_cap;
      }
    }

    if (lower_cap)
    {
      if (lower_cap_[*subcompid])
      {
        log(INFO, "WARNING: Multiple lower caps applied to %s::%s (component %d).  Taking maximum.",
                                  (*system_).name().c_str(), name().c_str(), *subcompid);
        *lower_cap_[*subcompid] = std::max(*lower_cap_[*subcompid], *lower_cap);
      }
      else
      {
        lower_cap_[*subcompid] = lower_cap;
      }
    }
     
  }
}

//*******************************************************************|************************************************************//
// fill the details of a bc assuming the buckettools schema and given a set of edge ids in a mesh function
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_bc_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string bcname;                                                // bc name
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), bcname); 
  spud_err(buffer.str(), serr);
  
  buffer.str(""); buffer << optionpath << "/sub_components";         // how many subcomponents exist?
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsubcomp; i++)                                // loop over the subcomponents
  {
    buffer.str(""); buffer << optionpath << "/sub_components[" 
                                                        << i << "]";
    initialize_bc_component_(buffer.str(), bcname);                  // initialize this component
  }
}

//*******************************************************************|************************************************************//
// initialize the expression for a component of a bc assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_bc_component_(const std::string &optionpath,
                                                  const std::string &bcname)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;                                      // a buffer to assemble this bc's component name

  std::string bctype;                                                // get the bc type
  buffer.str(""); buffer << optionpath << "/type/name";
  serr = Spud::get_option(buffer.str(), bctype);
  spud_err(buffer.str(), serr);

  if (bctype=="Dirichlet")
  {

    std::string compname;
    buffer.str(""); buffer << optionpath << "/name";
    serr = Spud::get_option(buffer.str(), compname);                 // get the name of these subcomponents
    spud_err(buffer.str(), serr);
    namebuffer.str(""); namebuffer << bcname << compname;            // assemble a name for the bc expression

    GenericFunction_ptr bcfunction = fetch_bcfunction(namebuffer.str());// fetch the bc expression in the function bucket maps
    Expression_ptr bcexp = std::dynamic_pointer_cast<dolfin::Expression>(bcfunction);
    buffer.str(""); buffer << optionpath << "/type::Dirichlet";      // initialize the expression describing the bc component
                                                                     // (assuming it is of type Dirichlet)
    initialize_expression_(bcexp, buffer.str(), 
                           namebuffer.str());

  }
  else
  {
    tf_err("Unknown bc type.", "bctype = ", bctype.c_str());
  }
       
}

//*******************************************************************|************************************************************//
// initialize a coefficient that is of type expression or constant
//*******************************************************************|************************************************************//
void SpudFunctionBucket::allocate_coeff_expression_()
{
  std::stringstream buffer;                                          // optionpath buffer

  if (type_=="Expression" || type_=="Constant")                      // only proceed for expressions and constants
  {

    buffer.str(""); buffer << optionpath() << "/type/rank/value";
    function_ = allocate_expression_over_regions_(buffer.str(), (*(*system_).bucket()).current_time_ptr());
    oldfunction_ = allocate_expression_over_regions_(buffer.str(), (*(*system_).bucket()).old_time_ptr());
    iteratedfunction_ = function_;

    buffer.str(""); buffer << (*system_).name() << "::" << name();     // rename the function as SystemName::CoefficientName
    (*function_).rename(buffer.str(), buffer.str());

    buffer.str(""); buffer << (*system_).name() << "::Old" << name();  // rename the old function as SystemName::OldCoefficientName
    (*oldfunction_).rename(buffer.str(), buffer.str());

    buffer.str(""); buffer << optionpath() << "/type/rank/value[0]/functional";
    if(Spud::have_option(buffer.str()))
    {
      fill_constantfunctional_();
    }

    fill_is_();

    outputfunctionspace_ = (*(*system()).bucket()).fetch_visfunctionspace((*system()).mesh());

  }

}

//*******************************************************************|************************************************************//
// fill the details of a constant functional assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::fill_constantfunctional_()
{
  constantfunctional_ = ufc_fetch_constant_functional((*system_).name(),// get a pointer to the functional form from the ufc
                        name(), (*system_).mesh());
  (*constantfunctional_).set_cell_domains((*system_).celldomains());
  (*constantfunctional_).set_interior_facet_domains((*system_).facetdomains());
  (*constantfunctional_).set_exterior_facet_domains((*system_).facetdomains());

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // functional because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
  uint ncoeff = (*constantfunctional_).num_coefficients();           // how many coefficients does this functional require?
  for (uint i = 0; i < ncoeff; i++)
  {
    std::string uflsymbol = 
                          (*constantfunctional_).coefficient_name(i);// what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
    if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))    // a base ufl symbol was only inserted into the parent bucket's
    {                                                                // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
      std::string baseuflsymbol =                                    // what is the base ufl symbol?
            (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol);   // have we already registered a functionspace for this base ufl
                                                                     // symbol?
      if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
      {                                                              // no...
        FunctionSpace_ptr coefficientspace;
        coefficientspace = 
                         ufc_fetch_coefficientspace_from_constant_functional( // take a pointer to the functionspace from the ufc
                                      (*system_).name(), name(), 
                                      baseuflsymbol, 
                                      (*system_).mesh(),
                                      (*system_).periodicmap(),
                                      (*system_).facetdomains(),
                                      (*system_).masterids(),
                                      (*system_).slaveids());
        (*(*system_).bucket()).register_coefficientspace(            // and register it in the parent bucket's map
                                      coefficientspace, 
                                      baseuflsymbol);
      }
    }

  }
}

//*******************************************************************|************************************************************//
// allocate an expression over region ids from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
Expression_ptr SpudFunctionBucket::allocate_expression_over_regions_(
                                       const std::string &optionpath,
                                       const double_ptr time)
{
  return allocate_expression_over_regions_(optionpath, time, NULL);
}

//*******************************************************************|************************************************************//
// allocate an expression over region ids from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
Expression_ptr SpudFunctionBucket::allocate_expression_over_regions_(
                                       const std::string &optionpath,
                                       const double_ptr time,
                                       bool *time_dependent)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  Expression_ptr expression;

  int nvs = Spud::option_count(optionpath);
  if (nvs > 1)
  {
    if (time_dependent)
    {
      *time_dependent = false;
    }

    std::map< std::size_t, Expression_ptr > expressions;

    int lrank;
    std::vector< std::size_t > shape;

    for (uint i = 0; i < nvs; i++)
    {
      std::string expressionname;
      buffer.str(""); buffer << optionpath << "[" << i << "]/name";
      serr = Spud::get_option(buffer.str(), expressionname);         // get the name of this expression
      
      buffer.str(""); buffer << optionpath << "[" << i << "]";
      Expression_ptr tmpexpression;
      if (time_dependent)
      {
        bool tmp_time_dependent;
        tmpexpression = allocate_expression_(
                            buffer.str(), expressionname, time, 
                                                &tmp_time_dependent);
        *time_dependent = *time_dependent || tmp_time_dependent;
      }
      else
      {
        tmpexpression = allocate_expression_(buffer.str(), 
                                              expressionname, time);
      }

      if (i==0)                                                      // record the rank and shape of the expression
      {
        lrank = (*tmpexpression).value_rank();
        for (uint d = 0; d < lrank; d++)
        {
          shape.push_back((*tmpexpression).value_dimension(d));
        }
      }
      else                                                           // check its the same as previous expressions
      {
        assert(lrank==(*tmpexpression).value_rank());
        for (uint d = 0; d < lrank; d++)
        {
          assert(shape[d]==(*tmpexpression).value_dimension(d));
        }
      }

      buffer.str(""); buffer << optionpath << "[" << i << "]/region_ids";
        
      std::vector<int> regionids;
      serr = Spud::get_option(buffer.str(), regionids);              // get the list of region ids
      spud_err(buffer.str(), serr);
    
      for (std::vector<int>::const_iterator id = regionids.begin();
                            id != regionids.end(); id++)
      {
        
        size_t_Expression_it e_it = expressions.find((std::size_t)*id);              // check if this component already exists
        if (e_it != expressions.end())
        {
          tf_err("Expression for region id defined multiple times in expression.", "Region id: %d, expression: %s.", 
                 *id, name().c_str());
        }
        else
        {
          expressions[*id] = tmpexpression;                          // if it doesn't, insert it into the map
        }

      }

    }

    MeshFunction_size_t_ptr cell_ids = (*system_).celldomains();
    if (lrank==0)
    {
      expression.reset( new RegionsExpression(expressions, 
                                                        cell_ids) );
    }
    else
    {
      expression.reset( new RegionsExpression(shape, expressions, 
                                                        cell_ids) );
    }

  }
  else
  {
    std::string expressionname;
    buffer.str(""); buffer << optionpath << "[0]/name";
    serr = Spud::get_option(buffer.str(), expressionname);           // get the name of this expression (probably WholeMesh)

    buffer.str(""); buffer << optionpath << "[0]";

    expression = allocate_expression_(buffer.str(), expressionname,
                                      time, time_dependent);         // allocate the function from the optionpath

  }

  return expression;

}

//*******************************************************************|************************************************************//
// initialize an expression from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
Expression_ptr SpudFunctionBucket::allocate_expression_(
                                       const std::string &optionpath,
                                       const std::string &expressionname,
                                       const double_ptr time)
{
  return allocate_expression_(optionpath, expressionname, time, NULL);
}

//*******************************************************************|************************************************************//
// initialize an expression from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
Expression_ptr SpudFunctionBucket::allocate_expression_(
                                       const std::string &optionpath,
                                       const std::string &expressionname,
                                       const double_ptr time,
                                       bool *time_dependent)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  Expression_ptr expression;                                         // declare the pointer that will be returned
  
  std::stringstream constbuffer, pybuffer, funcbuffer, cppbuffer, intbuffer;// some string assembly streams
  constbuffer.str(""); constbuffer << optionpath << "/constant";     // for a constant
  pybuffer.str("");    pybuffer    << optionpath << "/python";       // for a python function
  cppbuffer.str("");   cppbuffer   << optionpath << "/cpp";          // for a cpp expression
  intbuffer.str("");   intbuffer   << optionpath << "/internal";     // for an internal expression
  funcbuffer.str("");  funcbuffer  << optionpath << "/functional";   // for a functional
  
  if (Spud::have_option(pybuffer.str()))                             // python requested
  {
    std::string pyfunction;                                          // the python function string
    serr = Spud::get_option(pybuffer.str(), pyfunction); 
    spud_err(pybuffer.str(), serr);
    
                                                                     // rank of a python function isn't in the default spud base
                                                                     // language so have added it... but it comes out as a string 
                                                                     // of course!
    std::string rankstring;                                          // bit of a hack
    buffer.str(""); buffer << pybuffer.str() << "/rank";
    serr = Spud::get_option(buffer.str(), rankstring); 
    spud_err(buffer.str(), serr);
    
    int lrank;
    lrank = atoi(rankstring.c_str());
    if(lrank==0)                                                      // scalar
    {
      expression.reset(new PythonExpression(pyfunction, time));
    }
    else if (lrank==1)                                                // vector
    {
      expression.reset(new PythonExpression(size(), pyfunction, time));
    }
    else if (lrank==2)
    {
      expression.reset(new PythonExpression(shape_, pyfunction, time));
    }
    else
    {
      tf_err("Unknown rank for python expression in init_exp_.", "Rank: %d", lrank);
    }

    if (time_dependent)                                              // if we've asked if this expression is time dependent
    {                                                                // ... it may be
      *time_dependent = (*std::dynamic_pointer_cast< PythonExpression >(expression)).time_dependent();
    }

  }
  else if (Spud::have_option(funcbuffer.str()))                      // not much we can do for this case here
  {                                                                  // except check it's a scalar and declare a constant

                                                                     // rank of a functional isn't in the default spud base
                                                                     // language so have added it... but it comes out as a string 
                                                                     // of course!
    std::string rankstring;                                          // bit of a hack
    buffer.str(""); buffer << funcbuffer.str() << "/rank";
    serr = Spud::get_option(buffer.str(), rankstring); 
    spud_err(buffer.str(), serr);
    
    int lrank;
    lrank = atoi(rankstring.c_str());
    assert(lrank==0);

    expression.reset( new dolfin::Constant(0.0) );                   // can't set this to the correct value yet as it depends on
                                                                     // having all its coefficients attached so for now just set it
                                                                     // to zero

    if (time_dependent)                                              // if we've asked if this expression is time dependent
    {                                                                // ... assume it is (otherwise why would you be using a functional?)
      *time_dependent = true;
    }

  }
  else if (Spud::have_option(cppbuffer.str()))                       // not much we can do for this case here except call the constructor
  {

    std::string typestring;
    buffer.str(""); buffer << optionpath << "/type";                 // work out what type this expression is (initial_condition,
                                                                     // boundary_condition, value) based on an attribute included in
                                                                     // the schema.  This helps minimize the likelihood of namespace
                                                                     // overlap between different cpp expressions
    serr = Spud::get_option(buffer.str(), typestring); 
    spud_err(buffer.str(), serr);

    expression = cpp_fetch_expression((*system()).name(), name(), 
                                       typestring, expressionname, 
                                       size(), shape_, 
                                       (*system()).bucket(), 
                                       system(), time);

    if (time_dependent)                                              // if we've asked if this expression is time dependent
    {                                                                // ... check (but by default assume it is)
      buffer.str(""); buffer << optionpath << "cpp/time_independent";
      *time_dependent = !Spud::have_option(buffer.str());
    }

  }
  else if (Spud::have_option(intbuffer.str()))                       // not much we can do for this case here except call the constructor
  {

    std::string algoname;
    buffer.str(""); buffer << intbuffer.str() << "/algorithm/name";  // work out which internal algorithm we're using
    serr = Spud::get_option(buffer.str(), algoname);
    spud_err(buffer.str(), serr);

    if (algoname == "SemiLagrangian")
    {

      std::pair< std::string, std::pair< std::string, std::string > > function, velocity, outside;

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/lookup_function/system/name";
      serr = Spud::get_option(buffer.str(), function.first, (*system()).name());
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/lookup_function/field";
      if (Spud::have_option(buffer.str()))
      {
        buffer << "/name";
        function.second.first = "field";
        serr = Spud::get_option(buffer.str(), function.second.second);
        spud_err(buffer.str(), serr);
      }
      else
      {
        buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/lookup_function/coefficient/name";
        function.second.first = "coefficient";
        serr = Spud::get_option(buffer.str(), function.second.second);
        spud_err(buffer.str(), serr);
      }

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/velocity/system/name";
      serr = Spud::get_option(buffer.str(), velocity.first, (*system()).name());
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/velocity/field";
      if (Spud::have_option(buffer.str()))
      {
        buffer << "/name";
        velocity.second.first = "field";
        serr = Spud::get_option(buffer.str(), velocity.second.second);
        spud_err(buffer.str(), serr);
      }
      else
      {
        buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/velocity/coefficient/name";
        velocity.second.first = "coefficient";
        serr = Spud::get_option(buffer.str(), velocity.second.second);
        spud_err(buffer.str(), serr);
      }

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/outside_value/system/name";
      serr = Spud::get_option(buffer.str(), outside.first, (*system()).name());
      spud_err(buffer.str(), serr);

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/outside_value/field";
      if (Spud::have_option(buffer.str()))
      {
        buffer << "/name";
        outside.second.first = "field";
        serr = Spud::get_option(buffer.str(), velocity.second.second);
        spud_err(buffer.str(), serr);
      }
      else
      {
        buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/outside_value/coefficient/name";
        outside.second.first = "coefficient";
        serr = Spud::get_option(buffer.str(), outside.second.second);
        spud_err(buffer.str(), serr);
      }

                                                                     // rank of an internal function isn't in the default spud base
                                                                     // language so have added it... but it comes out as a string 
                                                                     // of course!
      std::string rankstring;                                        // bit of a hack
      buffer.str(""); buffer << intbuffer.str() << "/rank";
      serr = Spud::get_option(buffer.str(), rankstring); 
      spud_err(buffer.str(), serr);

      
      int lrank;
      lrank = atoi(rankstring.c_str());
      if (lrank==0)                                                   // scalar
      {
        expression.reset(new SemiLagrangianExpression(
                                               (*system()).bucket(), time, 
                                               function, velocity, outside));
      }
      else if (lrank==1)                                              // vector
      {
        expression.reset(new SemiLagrangianExpression(size(),
                                               (*system()).bucket(), time, 
                                               function, velocity, outside));
      }
      else if (lrank==2)                                              // vector
      {
        expression.reset(new SemiLagrangianExpression(shape_,
                                               (*system()).bucket(), time, 
                                               function, velocity, outside));
      }
      else
      {
        tf_err("Unknown rank for semi lagrangian expression in init_exp_.", "Rank: %d", lrank);
      }

      if (time_dependent)                                            // if we've asked if this expression is time dependent
      {                                                              // ... assume it is (otherwise why would you be using an algorithm?)
        *time_dependent = true;
      }

    }
    else if (algoname != "Reference")
    {
      tf_err("Unknown algorithm name in allocate_expression_.", "Algorithm: %s", algoname.c_str());
    }

  }
  else if (Spud::have_option(constbuffer.str()))                     // finally the constant case
  {
    int lrank;
    serr = Spud::get_option_rank(constbuffer.str(), lrank);          // find out the rank in the schema
    spud_err(constbuffer.str(), serr);
    
    if(lrank==0)                                                     // scalar
    {
      double value;
      serr = Spud::get_option(constbuffer.str(), value);             // take it from the options
      spud_err(constbuffer.str(), serr);
      expression.reset( new dolfin::Constant(value) );
    }
    else if (lrank==1)                                               // vector
    {
      std::vector<double> values;
      serr = Spud::get_option(constbuffer.str(), values); 
      spud_err(constbuffer.str(), serr);
      assert(values.size()==size());
      expression.reset(new dolfin::Constant(values));
    }
    else if (lrank==2)
    {
      std::vector< std::vector<double> > values_arr; 
      std::vector<int> value_shape_int;
      serr = Spud::get_option_shape(constbuffer.str(), value_shape_int); spud_err(constbuffer.str(), serr);
      serr = Spud::get_option(constbuffer.str(), values_arr); spud_err(constbuffer.str(), serr);

      std::vector<std::size_t> value_shape(2);
      std::vector<double> values;
      value_shape[0] = value_shape_int[0];
      value_shape[1] = value_shape_int[1];
      for (std::vector< std::vector<double> >::const_iterator val = values_arr.begin(); val != values_arr.end(); val++)
      {
        values.insert(values.end(), (*val).begin(), (*val).end());
      }
      assert(values.size()==size());

      expression.reset(new dolfin::Constant(value_shape, values));
    }
    else
    {
      tf_err("Unknown rank for constant in init_exp_.", "Rank: %d", lrank);
    }

    if (time_dependent)                                              // if we've asked if this expression is time dependent
    {                                                                // ... it's not
      *time_dependent = false;
    }

  } 
  
  return expression;                                                 // return the allocated expression
  
}

//*******************************************************************|************************************************************//
// take a reference to a function from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
GenericFunction_ptr SpudFunctionBucket::take_function_reference_(
                                       const std::string &optionpath,
                                       const std::string &expressionname,
                                       const double_ptr time)
{
  return take_function_reference_(optionpath, expressionname, time, NULL);
}

//*******************************************************************|************************************************************//
// take a reference to a function from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
GenericFunction_ptr SpudFunctionBucket::take_function_reference_(
                                       const std::string &optionpath,
                                       const std::string &expressionname,
                                       const double_ptr time,
                                       bool *time_dependent)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  GenericFunction_ptr function;                                      // declare the pointer that will be returned
  
  std::stringstream intbuffer;                                       // some string assembly streams
  intbuffer.str("");   intbuffer   << optionpath << "/internal";     // for an internal expression
  
  if (Spud::have_option(intbuffer.str()))
  {

    std::string algoname;
    buffer.str(""); buffer << intbuffer.str() << "/algorithm/name";  // work out which internal algorithm we're using
    serr = Spud::get_option(buffer.str(), algoname);
    spud_err(buffer.str(), serr);

    if (algoname == "Reference")
    {

      std::string systemname, functionname;

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/system/name";
      serr = Spud::get_option(buffer.str(), systemname, (*system()).name());
      spud_err(buffer.str(), serr);

      SystemBucket_ptr functionsystem = (*(*system()).bucket()).fetch_system(systemname);
      FunctionBucket_ptr functionbucket;

      buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/field";
      if (Spud::have_option(buffer.str()))
      {
        buffer << "/name";
        serr = Spud::get_option(buffer.str(), functionname);
        spud_err(buffer.str(), serr);

        functionbucket = (*functionsystem).fetch_field(functionname);
      }
      else
      {
        buffer.str(""); buffer << intbuffer.str() 
                        << "/algorithm/coefficient/name";
        serr = Spud::get_option(buffer.str(), functionname);
        spud_err(buffer.str(), serr);

        functionbucket = (*functionsystem).fetch_coeff(functionname);
      }

                                                                     // rank of an internal function isn't in the default spud base
                                                                     // language so have added it... but it comes out as a string 
                                                                     // of course!
      std::string rankstring;                                        // bit of a hack
      buffer.str(""); buffer << intbuffer.str() << "/rank";
      serr = Spud::get_option(buffer.str(), rankstring); 
      spud_err(buffer.str(), serr);

      
      int lrank;
      lrank = atoi(rankstring.c_str());

      const std::size_t functionrank = (*functionbucket).rank();

      if (functionrank != lrank)
      {
        tf_err("Rank of referenced function does not match rank of referee.", 
               "Reference rank: %d, Referee rank: %d", 
               functionrank, lrank);
      }

      function = (*functionbucket).genericfunction_ptr(time);

      if (time_dependent)                                            // clearly this could be a time-dependent function but it is
      {                                                              // not the responsibility of this functionbucket to update
        *time_dependent = false;                                     // it if it is so we record it as not being time-dependent!
      }

    }
  }
  
  return function;                                                    // return the function reference
  
}

//*******************************************************************|************************************************************//
// allocate an expression over region ids from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_expression_over_regions_(
                                       Expression_ptr expression,
                                       const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  int nvs = Spud::option_count(optionpath);
  if (nvs > 1)
  {
    std::map< std::size_t, Expression_ptr > expressions = 
      (*std::dynamic_pointer_cast< RegionsExpression >(expression)).expressions();

    for (uint i = 0; i < nvs; i++)
    {
      std::string expressionname;
      buffer.str(""); buffer << optionpath << "[" << i << "]/name";
      serr = Spud::get_option(buffer.str(), expressionname);         // get the name of this expression
      
      buffer.str(""); buffer << optionpath << "[" << i << "]/region_ids";
      std::vector<int> regionids;
      serr = Spud::get_option(buffer.str(), regionids);              // get the list of region ids
      spud_err(buffer.str(), serr);
    
      size_t_Expression_const_it e_it = expressions.find(regionids[0]);
      if (e_it == expressions.end())
      {
        tf_err("Unknown region id in RegionsExpression map.", "Region id: %d", regionids[0]);
      }
      else
      {
        buffer.str(""); buffer << optionpath << "[" << i << "]";
        initialize_expression_((*e_it).second, buffer.str(), expressionname);
      }

    }

  }
  else
  {
    std::string expressionname;
    buffer.str(""); buffer << optionpath << "[0]/name";
    serr = Spud::get_option(buffer.str(), expressionname);           // get the name of this expression (probably WholeMesh)

    buffer.str(""); buffer << optionpath << "[0]";

    initialize_expression_(expression, buffer.str(), expressionname);// initialize the expression from the optionpath

  }

}

//*******************************************************************|************************************************************//
// initialize an expression from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_expression_(
                                       Expression_ptr expression,
                                       const std::string &optionpath,
                                       const std::string &expressionname)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error

  std::stringstream funcbuffer, cppbuffer, intbuffer;                // some string assembly streams
  cppbuffer.str("");   cppbuffer   << optionpath << "/cpp";          // for a cpp expression
  intbuffer.str("");   intbuffer   << optionpath << "/internal";     // for an internal expression
  funcbuffer.str("");  funcbuffer  << optionpath << "/functional";   // for a functional

                                                                     // nothing to do here for python or constant as they are
                                                                     // already automaticall initialized... magic!

  if (Spud::have_option(funcbuffer.str()))                           // constant functionals
  {                                      

    assert(constantfunctional_);                                     // this should be valid as this option is only available for
                                                                     // Constant coefficients (i.e. not with multiple regions or
                                                                     // for boundary or initial conditions), hence we should always
                                                                     // be evaluating a value expression...
    std::string typestring;
    buffer.str(""); buffer << optionpath << "/type";                 // work out what type this expression is for assertion
    serr = Spud::get_option(buffer.str(), typestring); 
    spud_err(buffer.str(), serr);

    assert(typestring=="value");                                     // double check it's a value expression

    double value = dolfin::assemble(*constantfunctional_);
    *std::dynamic_pointer_cast< dolfin::Constant >(function_) = value;

  }
  else if (Spud::have_option(cppbuffer.str()))                       // cpp expressions
  {

    std::string typestring;
    buffer.str(""); buffer << optionpath << "/type";                 // work out what type this expression is (initial_condition,
                                                                     // boundary_condition, value) based on an attribute included in
                                                                     // the schema.  This helps minimize the likelihood of namespace
                                                                     // overlap between different cpp expressions
    serr = Spud::get_option(buffer.str(), typestring); 
    spud_err(buffer.str(), serr);

    cpp_init_expression(expression, (*system()).name(), name(), 
                                     typestring, expressionname);

  }
  else if (Spud::have_option(intbuffer.str()))                       // not much we can do for this case here except call the constructor
  {

    std::string algoname;
    buffer.str(""); buffer << intbuffer.str() << "/algorithm/name";  // work out which internal algorithm we're using
    serr = Spud::get_option(buffer.str(), algoname);
    spud_err(buffer.str(), serr);

    if (algoname == "SemiLagrangian")
    {
      (*std::dynamic_pointer_cast< SemiLagrangianExpression >(expression)).init();
    }
    else if (algoname != "Reference")
    {
      tf_err("Unknown algorithm name in initialize_expression_.", "Algorithm: %s", algoname.c_str());
    }

  }
  
}

//*******************************************************************|************************************************************//
// checkpoint the options file
//*******************************************************************|************************************************************//
void SpudFunctionBucket::checkpoint_options_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;

  namebuffer.str(""); namebuffer << (*(*system()).bucket()).output_basename() << "_" 
                                 << (*system()).name() << "_" 
                                 << (*(*system()).bucket()).checkpoint_count();
  buffer.str(""); buffer << optionpath()
                                  << "/type[0]/rank[0]/initial_condition";
  int nics = Spud::option_count(buffer.str());
  for (int i = nics-1; i >= 0; i--)
  {
    buffer.str(""); buffer << optionpath() << "/type[0]/rank[0]/initial_condition[" << i << "]";
    serr = Spud::delete_option(buffer.str());
    spud_err(buffer.str(), serr);
  }

  buffer.str(""); buffer << optionpath()
                                  << "/type[0]/rank[0]/initial_condition::WholeMesh/file";
  serr = Spud::set_option(buffer.str(), namebuffer.str());
  spud_err_accept(buffer.str(), serr, Spud::SPUD_NEW_KEY_WARNING);

  buffer.str(""); buffer << optionpath()
                                  << "/type[0]/rank[0]/initial_condition::WholeMesh/type";
  serr = Spud::set_option_attribute(buffer.str(), "initial_condition");
  spud_err_accept(buffer.str(), serr, Spud::SPUD_NEW_KEY_WARNING);

  buffer.str(""); buffer << optionpath()
                                  << "/type[0]/rank[0]/initial_condition::WholeMesh/file/__value/lines";
  serr = Spud::set_option_attribute(buffer.str(), "1");
  spud_err_accept(buffer.str(), serr, Spud::SPUD_NEW_KEY_WARNING);

}

