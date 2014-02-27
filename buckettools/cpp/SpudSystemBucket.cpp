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


#include "SpudSystemBucket.h"
#include "BoostTypes.h"
#include "SystemSolversWrapper.h"
#include "SpudBase.h"
#include "SpudFunctionBucket.h"
#include "SpudSolverBucket.h"
#include "PythonPeriodicMap.h"
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

  fill_periodicbcs_();                                              // fill in data about periodic bcs on the system

  fill_systemfunction_();                                            // register the functionspace and system functions

  fill_fields_();                                                    // initialize the fields (subfunctions) of this system

  fill_dirichletbcs_();                                              // fill in data about the bcs relative to the system (excludes
                                                                     // periodic)
  fill_points_();                                                    // initialize the reference point array of this system
 
  fill_coeffs_();                                                    // initialize the coefficient expressions (and constants)
                                                                     // (can't do coefficient functions now because it's unlikely we 
                                                                     // have all the coefficient functionspaces)

  fill_solvers_();                                                   // initialize the nonlinear solvers in this system

}

//*******************************************************************|************************************************************//
// loop over the coefficients and allocate any that are coefficient functions
//*******************************************************************|************************************************************//
void SpudSystemBucket::allocate_coeff_function()
{
  
  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();// loop over all coefficients
                                                              f_it++) 
  {                                                                  // recast as a spud derived class and initialize
    (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).allocate_coeff_function();
  }                                                                  // (check that this is a coefficient function within this
                                                                     // function)

}

//*******************************************************************|************************************************************//
// attach coefficients to forms and functionals
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_forms()
{
  dolfin::info("Attaching coeffs for system %s", name().c_str());
  attach_all_coeffs_();                                              // attach the coefficients to form and functionals
                                                                     // this happens here as some coefficients depend on functionals
                                                                     // to be evaluated
}

//*******************************************************************|************************************************************//
// attach coefficients to forms and functionals
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_fields_and_coefficients()
{
  dolfin::info("Initializing fields and coefficients for system %s", name().c_str());

  for (FunctionBucket_it f_it = fields_begin(); f_it != fields_end();
                                                              f_it++)
  {
    (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_field();
  } 

  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();
                                                              f_it++)
  {
    (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_coeff_expression();
  } 

                                                                     // after this point, we are allowed to start calling evals on
                                                                     // some of the expressions that have just been initialized
                                                                     // (this wasn't allowed up until now as all cpp expressions
                                                                     // potentially need initializing before eval will return the
                                                                     // right answer)
                                                                     // NOTE: even now there are potential inter dependencies! We
                                                                     // just deal with them by evaluating things in the order the
                                                                     // user specified.

  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();
                                                              f_it++)
  {
    (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_coeff_function();
  } 

  if (fields_size()>0)
  {
    apply_ic_();                                                     // apply the initial condition to the system function
    apply_dirichletbc_();                                            // apply the Dirichlet boundary conditions we just collected
    apply_referencepoints_();                                        // apply the reference points we just collected
  }

}

//*******************************************************************|************************************************************//
// initialize matrices described by this system's forms
//*******************************************************************|************************************************************//
void SpudSystemBucket::initialize_solvers()
{
  if (solve_location()!=SOLVE_NEVER)
  {
    dolfin::info("Initializing matrices for system %s", name().c_str());

    for (SolverBucket_it s_it = solvers_begin(); s_it != solvers_end();// loop over the solver buckets
                                                                s_it++)
    {
      (*boost::dynamic_pointer_cast< SpudSolverBucket >((*s_it).second)).initialize();
                                                                       // perform a preassembly of all the matrices to set up
                                                                       // sparsities etc. and set up petsc objects
    }
  }
}

//*******************************************************************|************************************************************//
// make a partial copy of the provided system bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SpudSystemBucket::copy_diagnostics(SystemBucket_ptr &system, Bucket_ptr &bucket) const
{

  if(!system)
  {
    system.reset( new SpudSystemBucket(optionpath_, &(*bucket)) );
  }

  SystemBucket::copy_diagnostics(system, bucket);

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
// checkpoint the options file
//*******************************************************************|************************************************************//
void SpudSystemBucket::checkpoint_options_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  if (solve_location()==SOLVE_START)                                 // do not solve again if this system was only meant to be
  {                                                                  // solved once
    std::string location = "never";
    buffer.str(""); buffer << optionpath() << "/solve/name";
    serr = Spud::set_option_attribute(buffer.str(), location);
    spud_err(buffer.str(), serr);
  }

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

  std::string location;
  buffer.str(""); buffer << optionpath() << "/solve/name";
  serr = Spud::get_option(buffer.str(), location);
  spud_err(buffer.str(), serr);
  if (location=="in_timeloop")
  {
    solve_location_ = SOLVE_TIMELOOP;
  }
  else if (location=="at_start")
  {
    solve_location_ = SOLVE_START;
  }
  else if (location=="with_diagnostics")
  {
    solve_location_ = SOLVE_DIAGNOSTICS;
  }
  else if (location=="never")
  {
    solve_location_ = SOLVE_NEVER;
  }
  else
  {
    dolfin::error("Unknown solve location for system %s.", name().c_str());
  }

  change_calculated_.reset( new bool(false) );                       // assume the change hasn't been calculated yet

  solved_.reset( new bool(false) );                                  // assume the system hasn't been solved yet

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

  functionspace_ = ufc_fetch_functionspace(name(), mesh(),           // fetch the first functionspace we can grab from the ufc 
                                           periodicmap(),            // for this system
                                           facetdomains(), 
                                           masterids(), slaveids());   
                                                                     

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
      size_t_Expression_it e_it = icexpressions.find(component);       // check if this component already exists
      if (e_it != icexpressions.end())
      {
        dolfin::error(                                               // if it does, issue an error
        "IC Expression with component number %d already exists in icexpressions map.", 
                                                          component);
      }
      else
      {
        icexpressions[component] = (*field).icexpression();          // if it doesn't, insert it into the map
      }

      component += (*(*field).icexpression()).value_size();          // increment the component count by the size of this field
                                                                     // (i.e. no. of scalar components)
    }
  }

  if (!icexpressions.empty())
  {
    collect_ics_(component, icexpressions);        // collect all the ics together into a new initial condition expression
  }

}

//*******************************************************************|************************************************************//
// fill in the data about the system bcs
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_periodicbcs_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str("");  buffer << optionpath() << "/boundary_condition";  // find out how many system bcs we have
  int nbcs = Spud::option_count(buffer.str());
  for (uint i = 0; i < nbcs; i++)                                    // loop over the system bcs in the options dictionary
  {

    std::vector<int> bcids;
    buffer.str(""); buffer << optionpath() << 
                    "/boundary_condition[" << i << 
                    "]/boundary_ids";
    serr = Spud::get_option(buffer.str(), bcids);
    spud_err(buffer.str(), serr);

    std::string bctype;
    buffer.str(""); buffer << optionpath() << 
                    "/boundary_condition[" << i << 
                    "]/sub_components[0]/type/name";
    serr = Spud::get_option(buffer.str(), bctype);
    spud_err(buffer.str(), serr);

    if(bctype=="Periodic")
    {
      for (uint j = 0; j < bcids.size(); j++)
      {
        masterids_.push_back(bcids[j]);
      }

      std::vector<int> slaveids_int;
      buffer.str(""); buffer << optionpath() << 
                      "/boundary_condition[" << i << 
                      "]/sub_components::All/type::Periodic/slave_boundary_ids";
      serr = Spud::get_option(buffer.str(), slaveids_int);
      spud_err(buffer.str(), serr);
      for (uint j = 0; j < slaveids_int.size(); j++)
      {
        slaveids_.push_back(slaveids_int[j]);
      }

      std::string function;
      buffer.str(""); buffer << optionpath() << 
                      "/boundary_condition[" << i << 
                      "]/sub_components::All/type::Periodic/coordinate_map";
      serr = Spud::get_option(buffer.str(), function);
      spud_err(buffer.str(), serr);

      periodicmap_.reset( new PythonPeriodicMap(function) );

    }
    else
    {
      dolfin::error("Unknown bc type.");
    }

  }

  buffer.str("");  buffer << "/io/debugging/periodic_boundaries";  // output debugging info?
  if (Spud::have_option(buffer.str()) && periodicmap())
  {
    buffer.str("");  buffer << (*bucket()).output_basename() << "_debugging_" << name() << "_periodic_boundaries.pvd";
    dolfin::File file(buffer.str());
    dolfin::MeshFunction<std::size_t> master_slave_entities;
    for (uint i = 0; i <= (*mesh()).topology().dim(); i++)
    {
      master_slave_entities = dolfin::PeriodicBoundaryComputation::masters_slaves(mesh(),
                                                      *periodicmap(), *facetdomains(), masterids(), slaveids(), i);
      file << master_slave_entities;
    }
  }
}

//*******************************************************************|************************************************************//
// fill in the data about the system bcs
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_dirichletbcs_()
{
  for (int_FunctionBucket_const_it f_it = orderedfields_begin();     // loop over all the fields
                                f_it != orderedfields_end(); f_it++)
  {
    for (int_DirichletBC_const_it                                    // loop over the dirichlet bcs
          b_it = (*(*f_it).second).ordereddirichletbcs_begin(); 
          b_it != (*(*f_it).second).ordereddirichletbcs_end(); b_it++)
    {
      dirichletbcs_.push_back(&(*boost::dynamic_pointer_cast<dolfin::DirichletBC>((*b_it).second)));                       // add the bcs to a std vector
    }
  }
}

//*******************************************************************|************************************************************//
// fill in the data about the system points (just grabs them from the fields)
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill_points_()
{
  for (int_FunctionBucket_const_it f_it = orderedfields_begin();     // loop over all the fields
                                f_it != orderedfields_end(); f_it++)
  {
    for (ReferencePoints_const_it                                     // loop over all the points
          p_it = (*(*f_it).second).points_begin(); 
          p_it != (*(*f_it).second).points_end(); p_it++)
    {
      points_.push_back((*p_it).second);                             // add the point to a std vector
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

