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


#include "Bucket.h"
#include "BoostTypes.h"
#include "SpudSystemBucket.h"
#include "SystemSolversWrapper.h"
#include "SpudBase.h"
#include "SpudFunctionBucket.h"
#include "SpudSolverBucket.h"
#include "SpudFunctionalBucket.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudSystemBucket::SpudSystemBucket(const std::string &optionpath, 
                                            Bucket* bucket) : 
                                            optionpath_(optionpath), 
                                            SystemBucket(bucket)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudSystemBucket::~SpudSystemBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// fill the system bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill()
{
  fill_base_();                                                      // fill in the base data (could be called from the
                                                                     // constructor?)

  fill_systemfunction_();                                            // register the functionspace and system functions

  fill_fields_();                                                    // initialize the fields (subfunctions) of this system

  fill_coeffs_();                                                    // initialize the coefficient expressions (and constants)
                                                                     // (can't do coefficient functions now because it's unlikely we 
                                                                     // have all the coefficient functionspaces)

  fill_solvers_();                                                   // initialize the nonlinear solvers in this system

  fill_functionals_();

}

//*******************************************************************|************************************************************//
// loop over the coefficients and allocate any that are coefficient functions
//*******************************************************************|************************************************************//
void SpudSystemBucket::allocate_coeff_function()
{
  
  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();// loop over all coefficients
                                                              f_it++) 
  {                                                                  // recast as a spud derived class and initialize
    (*(std::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).allocate_coeff_function();
  }                                                                  // (check that this is a coefficient function within this
                                                                     // function)

}

//*******************************************************************|************************************************************//
// allocate bcs
//*******************************************************************|************************************************************//
void SpudSystemBucket::allocate_bcs()
{

  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end();
                                                              f_it++)
  {
    (*(std::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).allocate_bcs();
  } 

  fill_bcs_();                                                       // fill in data about the bcs relative to the system (includes
                                                                     // periodic)

}

//*******************************************************************|************************************************************//
// attach coefficients to forms and functionals
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_fields_and_coefficient_expressions()
{
  log(INFO, "Initializing fields and coefficient expressions for system %s", name().c_str());

  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end();
                                                              f_it++)
  {
    (*(std::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_field();
  } 

  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();
                                                              f_it++)
  {
    (*(std::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_coeff_expression();
  } 

}

//*******************************************************************|************************************************************//
// attach coefficients to forms and functionals
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_coefficient_functions()
{
  log(INFO, "Initializing coefficient functions for system %s", name().c_str());

  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();
                                                              f_it++)
  {
    (*(std::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_coeff_function();
  } 

}

//*******************************************************************|************************************************************//
// initialize matrices described by this system's forms
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_solvers()
{
  log(INFO, "Initializing matrices for system %s", name().c_str());

  for (SolverBucket_it s_it = solvers_begin(); s_it != solvers_end();// loop over the solver buckets
                                                              s_it++)
  {
    if ((*(*s_it).second).solve_location() != SOLVE_NEVER)
    {
      (*std::dynamic_pointer_cast< SpudSolverBucket >((*s_it).second)).initialize();
                                                                     // perform a preassembly of all the matrices to set up
                                                                     // sparsities etc. and set up petsc objects
    }
  }
}

//*******************************************************************|************************************************************//
// return a vector of GenericFunctions that are to be included in the visualization file(s) from this system
//*******************************************************************|************************************************************//
std::vector< GenericFunction_ptr > SpudSystemBucket::collect_vis_functions() const
{
  std::vector< GenericFunction_ptr > functions;

  for (FunctionBucket_const_it f_it = fields_begin(); f_it != fields_end();
                                                              f_it++)
  {
    if ((*(*f_it).second).include_in_visualization())
    {
      functions.push_back( (*(*f_it).second).function() );
    }
    if ((*(*f_it).second).include_residual_in_visualization())
    {
      functions.push_back( (*(*f_it).second).residualfunction() );
    }
  }
 
  for (FunctionBucket_const_it c_it = coeffs_begin(); c_it != coeffs_end(); c_it++)
  {
    if ((*(*c_it).second).include_in_visualization())
    {
      functions.push_back( (*(*c_it).second).function() );
    }
  }

  return functions;
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the spud system
//*******************************************************************|************************************************************//
const std::string SpudSystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << " (" 
                                << optionpath() << ")" << std::endl;
  indent++;
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// fill the system bucket base data 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_base_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // get the system name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/ufl_symbol";           // get the system ufl symbol
  serr = Spud::get_option(buffer.str(), uflsymbol_); 
  spud_err(buffer.str(), serr);

  std::string meshname;                                              // get the system mesh
  buffer.str(""); buffer << optionpath() << "/mesh/name";
  serr = Spud::get_option(buffer.str(), meshname); 
  spud_err(buffer.str(), serr);
  mesh_ = (*bucket_).fetch_mesh(meshname);                           // and extract it from the bucket
  celldomains_ = (*bucket_).fetch_celldomains(meshname);             // along with the cell domains
  facetdomains_ = (*bucket_).fetch_facetdomains(meshname);           // and facet domains

  change_calculated_.reset( new bool(false) );                       // assume the change hasn't been calculated yet

}

//*******************************************************************|************************************************************//
// fill the system functionspace and function data 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_systemfunction_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str("");  buffer << optionpath() << "/field";               // find out how many fields we have
  int nfields = Spud::option_count(buffer.str());
  if (nfields==0)                                                    // and don't do anything if there are no fields
  {
    return;
  }

  functionspace_ = ufc_fetch_functionspace(name(), mesh());          // fetch the first functionspace we can grab from the ufc 
                                                                     // for this system

  function_.reset( new dolfin::Function(functionspace_) );           // declare the function on this functionspace
  buffer.str(""); buffer << name() << "::Function";
  (*function_).rename( buffer.str(), buffer.str() );

  oldfunction_.reset( new dolfin::Function(functionspace_) );        // declare the old function on this functionspace
  buffer.str(""); buffer << name() << "::OldFunction";
  (*oldfunction_).rename( buffer.str(), buffer.str() );

  iteratedfunction_.reset( new dolfin::Function(functionspace_) );   // declare the iterated function on this functionspace
  buffer.str(""); buffer << name() << "::IteratedFunction";
  (*iteratedfunction_).rename( buffer.str(), buffer.str() );

  changefunction_.reset( new dolfin::Function(functionspace_) );     // declare the change in the function between timesteps
  buffer.str(""); buffer << name() << "::TimestepChange";
  (*changefunction_).rename( buffer.str(), buffer.str() );

  residualfunction_.reset( new dolfin::Function(functionspace_) );   // declare the residual of the system as a function
  buffer.str(""); buffer << name() << "::Residual";
  (*residualfunction_).rename( buffer.str(), buffer.str() );

  if ((Spud::option_count(optionpath()+"/nonlinear_solver/type::SNES/monitors/visualization")+
       Spud::option_count(optionpath()+"/nonlinear_solver/type::SNES/monitors/convergence_file"))>0)
  {
    snesupdatefunction_.reset( new dolfin::Function(functionspace_) );
    buffer.str(""); buffer << name() << "::SNESUpdateFunction";
    (*snesupdatefunction_).rename( buffer.str(), buffer.str() );
  }

}

//*******************************************************************|************************************************************//
// fill in the data about each functional of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_functionals_()
{
  std::stringstream buffer;

  buffer.str(""); buffer << optionpath() << "/functional";
  int nfunctionals = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfunctionals; i++)
  {
    buffer.str(""); buffer << optionpath() << "/functional[" << i << "]";

    SpudFunctionalBucket_ptr functional(new SpudFunctionalBucket( buffer.str(), this ));
    (*functional).fill();
    register_functional(functional, (*functional).name());

  }
}

//*******************************************************************|************************************************************//
// fill in the data about each field (or subfunction) of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_fields_()
{
  std::stringstream buffer;                                          // optionpath buffer

                                                                     // prepare the system initial condition expression:
  uint component = 0;                                               // initialize a counter for the scalar components of this
                                                                     // system
  std::map< std::size_t, Expression_ptr > icexpressions;             // set up a map from scalar component to initial condition expression

  buffer.str("");  buffer << optionpath() << "/field";               // find out how many fields we have
  int nfields = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfields; i++)                                 // loop over the fields in the options dictionary
  {
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]";

                                                                     // declare a new field function bucket assuming this system is
                                                                     // its parent
    SpudFunctionBucket_ptr field(new SpudFunctionBucket( buffer.str(), this ));
    (*field).fill_field(i);                                          // fill in this field (providing its index in the system)
    register_field(field, (*field).name());                          // register this field in the system bucket
                                  
    if ((*field).icexpression())
    {
                                                                     // insert the field's initial condition expression into a 
                                                                     // temporary system map:
      size_t_Expression_it e_it = icexpressions.find(component);     // check if this component already exists
      if (e_it != icexpressions.end())
      {
        tf_err("IC expression with given component number already exists in inexpressions map.", "Component: %d", component);
      }
      else
      {
        icexpressions[component] = (*field).icexpression();          // if it doesn't, insert it into the map
      }
    }

    component += (*field).size();                                    // increment the component count by the size of this field
                                                                     // (i.e. no. of scalar components)
  }

  if (!icexpressions.empty())
  {
    collect_ics_(component, icexpressions);                          // collect all the ics together into a new initial condition expression
  }

}

//*******************************************************************|************************************************************//
// fill in the data about the system bcs
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_bcs_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  for (FunctionBucket_const_it f_it = fields_begin();     // loop over all the fields
                                f_it != fields_end(); f_it++)
  {
    for (DirichletBC_const_it                                    // loop over the bcs
          b_it = (*(*f_it).second).dirichletbcs_begin(); 
          b_it != (*(*f_it).second).dirichletbcs_end(); b_it++)
    {
      bcs_.push_back((*b_it).second);                   // add the bcs to a std vector
    }
    for (ReferencePoint_const_it                                    // loop over all the points
          p_it = (*(*f_it).second).referencepoints_begin(); 
          p_it != (*(*f_it).second).referencepoints_end(); p_it++)
    {
      bcs_.push_back((*p_it).second);                    // add the point to a std vector
    }
  }

}

//*******************************************************************|************************************************************//
// fill in the data about each coefficient expression (or constant) of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_coeffs_()
{
  std::stringstream buffer;                                          // optionpath buffer
  
  buffer.str("");  buffer << optionpath() << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());                    // find out how many coefficients we have
  for (uint i = 0; i < ncoeffs; i++)                                 // and loop over them
  {
    buffer.str(""); buffer << optionpath() << "/coefficient[" 
                                                        << i << "]";

                                                                     // initialize a new function bucket for this coefficient
                                                                     // (regardless of type!) assuming this system is its parent
    SpudFunctionBucket_ptr coeff( new SpudFunctionBucket( buffer.str(), this ) );
    (*coeff).fill_coeff(i);                                          // fill the coefficient (this won't do much for coefficient
                                                                     // functions)
    register_coeff(coeff, (*coeff).name());                          // register this coefficient in the system

  }

}

//*******************************************************************|************************************************************//
// fill in the data about each solver of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_solvers_()
{
  std::stringstream buffer;                                          // optionpath buffer
  
  buffer.str("");  buffer << optionpath() << "/nonlinear_solver";
  int nsolvers = Spud::option_count(buffer.str());                   // find out how many nonlinear solvers there are
  for (uint i = 0; i < nsolvers; i++)                                // loop over them
  {
    buffer.str(""); buffer << optionpath() << "/nonlinear_solver[" 
                                                        << i << "]";

                                                                     // initialize a new solver bucket assuming this system is its
                                                                     // parent
    SpudSolverBucket_ptr solver( new SpudSolverBucket( buffer.str(), this ) );
    (*solver).fill();                                                // fill in the data about this solver bucket
    register_solver(solver, (*solver).name());                       // register the solver bucket in the system
  }
}

