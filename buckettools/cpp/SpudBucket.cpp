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


#include "GlobalPythonInstance.h"
#include "PythonDetectors.h"
#include "SpudBucket.h"
#include "SpudSystemBucket.h"
#include "SpudSystemsSolverBucket.h"
#include "SpudBase.h"
#include "BoostTypes.h"
#include "BucketDolfinBase.h"
#include "PointDetectors.h"
#include "StatisticsFile.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket() : Bucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket(const std::string &name) : Bucket(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudBucket::SpudBucket(const std::string &name, const std::string &optionpath) : 
                                optionpath_(optionpath), Bucket(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudBucket::~SpudBucket()
{
}

//*******************************************************************|************************************************************//
// fill the bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudBucket::fill()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  fill_globalparameters_();

  buffer.str(""); buffer << optionpath() << "/geometry/dimension";   // geometry dimension set in the bucket to pass it down to all
  serr = Spud::get_option(buffer.str(), dimension_);                 // systems (we assume this is the length of things that do
  spud_err(buffer.str(), serr);                                      // not have them independently specified)

  fill_timestepping_();                                              // fill in the timestepping options (if there are any)

  buffer.str(""); buffer << optionpath() << "/geometry/mesh";        // put the meshes into the bucket
  int nmeshes = Spud::option_count(buffer.str());
  for (uint i = 0; i<nmeshes; i++)                                   // loop over the meshes defined in the options file
  {
    buffer.str(""); buffer << optionpath() << "/geometry/mesh[" 
                                                        << i << "]";
    fill_meshes_(buffer.str());
  }
  
  fill_output_();                                                    // fill in the output options (if there are any)
                                                                     // has to happen after fill_meshes_ as it outputs the mesh

  buffer.str(""); buffer << optionpath() << "/system";
  int nsystems = Spud::option_count(buffer.str());
  for (uint i = 0; i<nsystems; i++)                                  // loop over the systems registering the base uflsymbols of any
  {                                                                  // coefficient functions contained in them
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    fill_baseuflsymbols_(buffer.str());
  }

  for (uint i = 0; i<nsystems; i++)                                  // loop over the systems *again*, this time filling all the
  {                                                                  // systems data into the bucket data structures
    buffer.str(""); buffer << optionpath() << "/system[" << i << "]";
    fill_systems_(buffer.str());
  }

  fill_systemssolvers_();                                            // fill in information about the systems solvers

  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *third* time, this time filling
                                  sys_it != systems_end(); sys_it++) // in the data for the coefficient functions contained within
  {                                                                  // them
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).allocate_coeff_function();
  }                                                                  // we couldn't do this before because we might not have had
                                                                     // the right functionspace available
  
  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *fourth* time, this time filling
                                  sys_it != systems_end(); sys_it++) // in the data for the bcs
  {                                                                  // 
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).allocate_bcs();
  }                                                                  // we couldn't do this before now because we allow reference
                                                                     // bcs so everything had to be allocated
  
  fill_uflsymbols_();                                                // now all the functions in the systems are complete we can 
                                                                     // register them in the bucket so it's easy to attach them
                                                                     // to the forms and functionals

  fill_adaptivetimestepping_();                                      // fill in the adaptive timestepping options (if there are any)
                                                                     // this can also be done now because all relevant pointers should be
                                                                     // allocated

  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *fifth* time, attaching the
                                  sys_it != systems_end(); sys_it++) // coefficients to the forms and functionals
  {
    (*(*sys_it).second).initialize_forms();
  }
  
  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *sixth* time, initializing
                                  sys_it != systems_end(); sys_it++) // the values of any expressions, functionals or functions
  {                                                                  // used by fields or coefficients
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).initialize_fields_and_coefficient_expressions();
  }
  
                                                                     // after this point, we are allowed to start calling evals on
                                                                     // some of the expressions that have just been initialized
                                                                     // (this wasn't allowed up until now as all cpp expressions
                                                                     // potentially need initializing before eval will return the
                                                                     // right answer)
                                                                     // NOTE: even now there are potential inter dependencies! We
                                                                     // just deal with them by evaluating coefficient functions 
                                                                     // in the order the user specified followed by fields last.

  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *seventh* time, initializing
                       sys_it != systems_end(); sys_it++)            // the values of any coefficient functions
  {                                                                  
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).initialize_coefficient_functions();
  }
  
  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *eighth* time, initializing
                       sys_it != systems_end(); sys_it++)            // the values of the fields
  {
    (*(*sys_it).second).evaluate_initial_fields();
  }
  
  for (SystemBucket_it sys_it = systems_begin();                     // loop over the systems for a *eighth* time, preassembling
                                  sys_it != systems_end(); sys_it++) // the matrices
  {
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*sys_it).second)).initialize_solvers();
  }
  
  fill_detectors_();                                                 // put the detectors in the bucket

  fill_diagnostics_();                                               // this should be called last because it initializes the
                                                                     // diagnostic files, which must use a complete bucket

  log(INFO, str().c_str());

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a dolfin mesh in the bucket (and spudbucket) data maps with a spud optionpath
//*******************************************************************|************************************************************//
void SpudBucket::register_mesh(Mesh_ptr mesh, 
                               const std::string &name, 
                               std::string optionpath,
                               MeshFunction_size_t_ptr celldomains, 
                               MeshFunction_size_t_ptr facetdomains)
{
  Mesh_hash_it m_it = meshes_.get<om_key_hash>().find(name);                                 // check if a mesh with this name already exists
  if (m_it != meshes_.get<om_key_hash>().end())
  {
    tf_err("Mesh already exists in spudbucket.", "Mesh name: %s", name.c_str());
  }
  else
  {
    meshes_.insert(om_item<const std::string,Mesh_ptr>(name,mesh));                                  // if not register it in the map
    mesh_optionpaths_.insert(om_item<const std::string,std::string>(name,optionpath));                            // also register its optionpath
    celldomains_.insert(om_item<const std::string,MeshFunction_size_t_ptr>(name,celldomains));
    facetdomains_.insert(om_item<const std::string,MeshFunction_size_t_ptr>(name,facetdomains));
  }
}

//*******************************************************************|************************************************************//
// return a string containing the mesh optionpath in the spudbucket data maps
//*******************************************************************|************************************************************//
std::string SpudBucket::fetch_mesh_optionpath(const std::string &name)
{
  string_hash_it s_it = mesh_optionpaths_.get<om_key_hash>().find(name);                     // check if a mesh with this name exists
  if (s_it == mesh_optionpaths_.get<om_key_hash>().end())
  {
    tf_err("Mesh does not exist in spudbucket.", "Mesh name: %s", name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the optionpath
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::mesh_optionpaths_begin()
{
  return mesh_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::mesh_optionpaths_begin() const
{
  return mesh_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::mesh_optionpaths_end()
{
  return mesh_optionpaths_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::mesh_optionpaths_end() const
{
  return mesh_optionpaths_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a detector in the bucket (and spudbucket) data maps with a spud optionpath
//*******************************************************************|************************************************************//
void SpudBucket::register_detector(GenericDetectors_ptr detector, 
                               const std::string &name, 
                               std::string optionpath)
{
  GenericDetectors_hash_it d_it = detectors_.get<om_key_hash>().find(name);                  // check if a detector set with this name already exists
  if (d_it != detectors_.get<om_key_hash>().end())
  {
    tf_err("Detector set already exists in spudbucket.", "Detector set name: %s", name.c_str());
  }
  else
  {
    detectors_.insert(om_item<const std::string, GenericDetectors_ptr>(name, detector));                          // if not register it in the map
    detector_optionpaths_.insert(om_item<const std::string, std::string>(name,optionpath));                        // also register its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string containing the detector optionpath in the spudbucket data maps
//*******************************************************************|************************************************************//
std::string SpudBucket::fetch_detector_optionpath(const std::string &name)
{
  string_hash_it s_it = detector_optionpaths_.get<om_key_hash>().find(name);                 // check if a mesh with this name exists
  if (s_it == detector_optionpaths_.get<om_key_hash>().end())
  {
    tf_err("Detector set does not exist in spudbucket.", "Detector set name: %s", name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the optionpath
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::detector_optionpaths_begin()
{
  return detector_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::detector_optionpaths_begin() const
{
  return detector_optionpaths_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the detector_optionpaths_ map
//*******************************************************************|************************************************************//
string_it SpudBucket::detector_optionpaths_end()
{
  return detector_optionpaths_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the mesh_optionpaths_ map
//*******************************************************************|************************************************************//
string_const_it SpudBucket::detector_optionpaths_end() const
{
  return detector_optionpaths_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the mesh_optionpaths_ data structure
//*******************************************************************|************************************************************//
const std::string SpudBucket::meshes_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = mesh_optionpaths_begin(); 
                            s_it != mesh_optionpaths_end(); s_it++ )
  {
    s << indentation << "Mesh " << (*s_it).first 
                    << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the detector_optionpaths_ data structure
//*******************************************************************|************************************************************//
const std::string SpudBucket::detectors_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = detector_optionpaths_begin(); 
                            s_it != detector_optionpaths_end(); s_it++ )
  {
    s << indentation << "Detector set " << (*s_it).first 
                    << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// read the global parameters set in the schema
//*******************************************************************|************************************************************//
void SpudBucket::fill_globalparameters_() const
{
  std::stringstream buffer;
  Spud::OptionError serr;

  buffer.str(""); buffer << "/global_parameters/python";
  if (Spud::have_option(buffer.str()))
  {
    std::string function;
    serr = Spud::get_option(buffer.str(), function);
    spud_err(buffer.str(), serr);
    (*GlobalPythonInstance::instance()).run(function);
  }

  buffer.str(""); buffer << "/global_parameters/dolfin/allow_extrapolation";
  if (Spud::have_option(buffer.str()))
  {
    dolfin::parameters["allow_extrapolation"] = true;
  }

  buffer.str(""); buffer << "/global_parameters/dolfin/ghost_mode";
  if (Spud::have_option(buffer.str()))
  {
    buffer << "/name";
    std::string ghost_mode;
    serr = Spud::get_option(buffer.str(), ghost_mode);
    spud_err(buffer.str(), serr);
    dolfin::parameters["ghost_mode"] = ghost_mode;
  }

}

//*******************************************************************|************************************************************//
// fill in any timestepping data (except adaptive stuff) or set up dummy values instead (zero essentially)
//*******************************************************************|************************************************************//
void SpudBucket::fill_timestepping_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  
  timestep_count_.reset( new int );
  *timestep_count_ = 0;                                              // the number of timesteps taken

  start_time_.reset( new double );
  buffer.str(""); buffer << "/timestepping/current_time";            // get the current time
  serr = Spud::get_option(buffer.str(), *start_time_, 0.0);          // may be non-zero at the start of simulations (e.g.
  spud_err(buffer.str(), serr);                                      // checkpoints) but assumed zero for steady simulations

  current_time_.reset( new double(start_time()) );                   // initialize the current (n+1) and
  old_time_.reset( new double(start_time()) );                       // old (n) times

  buffer.str(""); buffer << "/timestepping";
  if (Spud::have_option(buffer.str()))
  {
    buffer.str(""); 
    buffer << "/timestepping/timestep/coefficient::Timestep/ufl_symbol";
    serr = Spud::get_option(buffer.str(), timestep_.first);
    spud_err(buffer.str(), serr);

    buffer.str(""); 
    double timestep_value;
    buffer << "/timestepping/timestep/coefficient::Timestep/type::Constant/rank::Scalar/value::WholeMesh/constant";
    serr = Spud::get_option(buffer.str(), timestep_value);
    spud_err(buffer.str(), serr);
    timestep_.second.reset( new dolfin::Constant(timestep_value) );

    buffer.str(""); buffer << "/timestepping/steady_state";
    if (Spud::have_option(buffer.str()))
    {
      if ((Spud::option_count("/system/functional/include_in_steady_state")+
           Spud::option_count("/system/field/diagnostics/include_in_steady_state"))==0)
      {
        tf_err("Reqested a steady state check but selected no fields or functionals to include.", 
               "No fields or functionals in steady state check.");
      }

      steadystate_tol_.reset( new double );                          // get the steady state tolerance
      buffer.str(""); buffer << "/timestepping/steady_state/tolerance";
      serr = Spud::get_option(buffer.str(), *steadystate_tol_); 
      spud_err(buffer.str(), serr);
    }

    buffer.str(""); buffer << "/timestepping/finish_time";           // get the finish time (assumed zero for steady simulations)
    if (Spud::have_option(buffer.str()))
    {
      finish_time_.reset( new double );
      serr = Spud::get_option(buffer.str(), *finish_time_); 
      spud_err(buffer.str(), serr);
    }
    else
    {
      number_timesteps_.reset( new int );
      buffer.str(""); buffer << "/timestepping/number_timesteps";    // get the finish time (assumed zero for steady simulations)
      serr = Spud::get_option(buffer.str(), *number_timesteps_); 
      spud_err(buffer.str(), serr);
    }

    buffer.str(""); buffer << "/timestepping/walltime_limit";
    if (Spud::have_option(buffer.str()))
    {
      walltime_limit_.reset( new double );
      serr = Spud::get_option(buffer.str(), *walltime_limit_);
      spud_err(buffer.str(), serr);
    }

  }
  else
  {
    number_timesteps_.reset( new int(1) );

    timestep_.first = "";
    timestep_.second.reset( new dolfin::Constant(0.0) );

  }
  
  
}

//*******************************************************************|************************************************************//
// fill in any adaptive timestepping data 
//*******************************************************************|************************************************************//
void SpudBucket::fill_adaptivetimestepping_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  
  buffer.str(""); buffer << "/timestepping/timestep/adaptive";
  if (Spud::have_option(buffer.str()))
  {
    
    
    buffer.str(""); buffer <<
       "/timestepping/timestep/adaptive/constraint";               // constraints on which to adapt the timestep (required)
    int nconstraints = Spud::option_count(buffer.str());
    for (uint i = 0; i < nconstraints; i++)
    {
      buffer.str(""); buffer << 
         "/timestepping/timestep/adaptive/constraint[" << i << "]";
      
      std::string sysname;
      serr = Spud::get_option(buffer.str()+"/system/name", sysname);
      spud_err(buffer.str()+"/system/name", serr);

      std::string functionname;
      FunctionBucket_ptr function;
      if (Spud::have_option(buffer.str()+"/field/name"))
      {
        serr = Spud::get_option(buffer.str()+"/field/name", functionname);
        spud_err(buffer.str()+"/field/name", serr);
        function = (*fetch_system(sysname)).fetch_field(functionname);
      }
      else
      {
        serr = Spud::get_option(buffer.str()+"/coefficient/name", functionname);
        spud_err(buffer.str()+"/coefficient/name", serr);
        function = (*fetch_system(sysname)).fetch_coeff(functionname);
      }

      double maxvalue;
      serr = Spud::get_option(buffer.str()+"/requested_maximum_value", maxvalue);
      spud_err(buffer.str()+"/requested_maximum_value", serr);

      timestep_constraints_.push_back(std::make_pair( function, maxvalue ));

    }

    buffer.str(""); buffer <<
       "/timestepping/timestep/adaptive/adapt_period_in_timesteps";// timestep adapt period in timesteps
    if(Spud::have_option(buffer.str()))
    {
      timestepadapt_period_timesteps_.reset( new int );
      serr = Spud::get_option(buffer.str(), *timestepadapt_period_timesteps_);
      spud_err(buffer.str(), serr);
    }
    
    buffer.str(""); buffer << 
                "/timestepping/timestep/adaptive/adapt_period";    // timestep adapt period
    if(Spud::have_option(buffer.str()))
    {
      timestepadapt_period_.reset( new double );
      serr = Spud::get_option(buffer.str(), *timestepadapt_period_);
      spud_err(buffer.str(), serr);

      timestepadapt_time_.reset( new double(start_time()) );
    }

    buffer.str(""); buffer << 
            "/timestepping/timestep/adaptive/increase_tolerance";  // timestep adapt period
    if(Spud::have_option(buffer.str()))
    {
      timestep_increasetol_.reset( new double );
      serr = Spud::get_option(buffer.str(), *timestep_increasetol_);
      spud_err(buffer.str(), serr);
    }

  }
  
}

//*******************************************************************|************************************************************//
// fill in any output data
//*******************************************************************|************************************************************//
void SpudBucket::fill_output_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud option error
  
  buffer.str(""); buffer << "/io/output_base_name";                  // get the output base name
  serr = Spud::get_option(buffer.str(), output_basename_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << "/io/dump_periods/visualization_period"; // visualization period
  if(Spud::have_option(buffer.str()))
  {
    visualization_period_.reset( new double );
    serr = Spud::get_option(buffer.str(), *visualization_period_);
    spud_err(buffer.str(), serr);

    visualization_dumptime_.reset( new double(start_time()) );
  }

  buffer.str(""); buffer 
            << "/io/dump_periods/visualization_period_in_timesteps"; // visualization period in timesteps
  if(Spud::have_option(buffer.str()))
  {
    visualization_period_timesteps_.reset( new int );
    serr = Spud::get_option(buffer.str(), *visualization_period_timesteps_);
    spud_err(buffer.str(), serr);
  }

  visualization_count_.reset( new int(0) );
  
  buffer.str(""); buffer << "/io/dump_periods/statistics_period";    // statistics period
  if(Spud::have_option(buffer.str()))
  {
    statistics_period_.reset( new double );
    serr = Spud::get_option(buffer.str(), *statistics_period_);
    spud_err(buffer.str(), serr);

    statistics_dumptime_.reset( new double(start_time()) );
  }

  buffer.str(""); buffer 
            << "/io/dump_periods/statistics_period_in_timesteps";    // statistics period in timesteps
  if(Spud::have_option(buffer.str()))
  {
    statistics_period_timesteps_.reset( new int );
    serr = Spud::get_option(buffer.str(), *statistics_period_timesteps_);
    spud_err(buffer.str(), serr);
  }
  
  buffer.str(""); buffer << "/io/dump_periods/steady_state_period";  // steady state period
  if(Spud::have_option(buffer.str()))
  {
    steadystate_period_.reset( new double );
    serr = Spud::get_option(buffer.str(), *steadystate_period_);
    spud_err(buffer.str(), serr);

    steadystate_dumptime_.reset( new double(start_time()) );
  }

  buffer.str(""); buffer 
            << "/io/dump_periods/steady_state_period_in_timesteps";  // steady state period in timesteps
  if(Spud::have_option(buffer.str()))
  {
    steadystate_period_timesteps_.reset( new int );
    serr = Spud::get_option(buffer.str(), *steadystate_period_timesteps_);
    spud_err(buffer.str(), serr);
  }
  
  buffer.str(""); buffer << "/io/dump_periods/detectors_period";     // detectors period
  if(Spud::have_option(buffer.str()))
  {
    detectors_period_.reset( new double );
    serr = Spud::get_option(buffer.str(), *detectors_period_);
    spud_err(buffer.str(), serr);

    detectors_dumptime_.reset( new double(start_time()) );
  }

  buffer.str(""); buffer 
            << "/io/dump_periods/detectors_period_in_timesteps";     // detectors period in timesteps
  if(Spud::have_option(buffer.str()))
  {
    detectors_period_timesteps_.reset( new int );
    serr = Spud::get_option(buffer.str(), *detectors_period_timesteps_);
    spud_err(buffer.str(), serr);
  }
  
  buffer.str(""); buffer << "/io/checkpointing/checkpoint_period";   // checkpoint period
  if(Spud::have_option(buffer.str()))
  {
    checkpoint_period_.reset( new double );
    serr = Spud::get_option(buffer.str(), *checkpoint_period_);
    spud_err(buffer.str(), serr);

    checkpoint_time_.reset( new double(start_time()) );
    checkpoint_count_.reset( new int(0) );
  }

  buffer.str(""); buffer 
            << "/io/checkpointing/checkpoint_period_in_timesteps";   // checkpoint period in timesteps
  if(Spud::have_option(buffer.str()))
  {
    checkpoint_period_timesteps_.reset( new int );
    serr = Spud::get_option(buffer.str(), *checkpoint_period_timesteps_);
    spud_err(buffer.str(), serr);

    checkpoint_count_.reset( new int(0) );
  }

  if(Spud::option_count("/system/functional/output_cell_function") +
     Spud::option_count("/system/functional/output_facet_function") > 0)
  {
    int nmeshes = Spud::option_count("/geometry/mesh");
    for (uint i = 0; i<nmeshes; i++)                                   // loop over the meshes defined in the options file
    {
      buffer.str(""); buffer << optionpath() << "/geometry/mesh[" 
                                                          << i << "]/source/name";
      std::string source;                                                // get the source of the mesh (i.e. from file or internal)
      serr = Spud::get_option(buffer.str(), source); 
      spud_err(buffer.str(), serr);
      if(source!="File")
      {
        buffer.str(""); buffer << optionpath() << "/geometry/mesh[" 
                                                            << i << "]/name";
        std::string meshname;                                              // get the name of the mesh
        serr = Spud::get_option(buffer.str(), meshname); 
        spud_err(buffer.str(), serr);

        Mesh_ptr mesh = fetch_mesh(meshname);

        buffer.str(""); buffer << output_basename() << "_" << meshname << ".xml";
        dolfin::File meshfile(buffer.str());
        meshfile << *mesh;
      }
    }

  }

}

//*******************************************************************|************************************************************//
// create or get from a file a dolfin mesh object and insert it into the bucket data structures
//*******************************************************************|************************************************************//
void SpudBucket::fill_meshes_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string meshname;                                              // get the name of the mesh
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), meshname); 
  spud_err(buffer.str(), serr);
  
  std::string source;                                                // get the source of the mesh (i.e. from file or internal)
  buffer.str(""); buffer << optionpath << "/source/name";
  serr = Spud::get_option(buffer.str(), source); 
  spud_err(buffer.str(), serr);

  Mesh_ptr mesh;                                                     // initialize the pointer
  MeshFunction_size_t_ptr edgeids;
  MeshFunction_size_t_ptr cellids;

  if (source=="File")                                                // source is a file
  {
    std::string basename;                                            // get the base file name (without the .xml)
    buffer.str(""); buffer << optionpath << "/source/file";
    serr = Spud::get_option(buffer.str(), basename); 
    spud_err(buffer.str(), serr);


    std::ifstream file;                                              // dummy file stream to test if files exist
                                                                     // (better way of doing this?)
    
    std::stringstream filename;

    filename.str(""); filename << basename << ".xml";
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)
    {
      file.close();
      mesh.reset(new dolfin::Mesh(filename.str()));
    }
    else
    {
      filename.str(""); filename << basename << ".xml.gz";
      file.open(filename.str().c_str(), std::ifstream::in);
      if (file)
      {
        file.close();
        mesh.reset(new dolfin::Mesh(filename.str()));
      }
      else
      {
        tf_err("Could not find requested mesh.", 
               "%s.xml or %s.xml.gz not found.", basename.c_str(), basename.c_str());
      }
    }
    (*mesh).init();                                                  // initialize the mesh (maps between dimensions etc.)

    filename.str(""); filename << basename << "_facet_region.xml";   // check if the edge subdomain mesh function file exists
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)                                                        // if it does then attach it to the dolfin MeshData structure 
    {                                                                // using the dolfin reserved name for exterior facets
      file.close();

      edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, filename.str()));
    }
    else
    {
      filename.str(""); filename << basename << "_facet_region.xml.gz";// check if the edge subdomain mesh function file exists
      file.open(filename.str().c_str(), std::ifstream::in);
      if (file)                                                      // if it does then attach it to the dolfin MeshData structure 
      {                                                              // using the dolfin reserved name for exterior facets
        file.close();

        edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, filename.str()));
      }
    }

    filename.str(""); filename << basename << "_physical_region.xml";// check if the region subdomain mesh function file exists
    file.open(filename.str().c_str(), std::ifstream::in);
    if (file)                                                        // if it does then attach it to the dolfin MeshData structure 
    {                                                                // using the dolfin reserved name for cell domains
      file.close();

      cellids.reset(new dolfin::MeshFunction<std::size_t>(mesh, filename.str()));
    }
    else
    {
      filename.str(""); filename << basename << "_physical_region.xml.gz";// check if the region subdomain mesh function file exists
      file.open(filename.str().c_str(), std::ifstream::in);
      if (file)                                                        // if it does then attach it to the dolfin MeshData structure 
      {                                                                // using the dolfin reserved name for cell domains
        file.close();

        cellids.reset(new dolfin::MeshFunction<std::size_t>(mesh, filename.str()));
      }
    }

  }
  else if (source=="UnitInterval")                                   // source is an internally generated dolfin mesh
  {
    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitIntervalMesh(cells) );

    Side left(0, 0.0);
    Side right(0, 1.0);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 0, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
  }
  else if (source=="Interval")                                       // source is an internally generated dolfin mesh
  {
    double leftx;
    buffer.str(""); buffer << optionpath << "/source/left";
    serr = Spud::get_option(buffer.str(), leftx); 
    spud_err(buffer.str(), serr);

    double rightx;
    buffer.str(""); buffer << optionpath << "/source/right";
    serr = Spud::get_option(buffer.str(), rightx); 
    spud_err(buffer.str(), serr);

    int cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::IntervalMesh(cells, leftx, rightx) );

    Side left(0, leftx);
    Side right(0, rightx);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 0, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
  }
  else if (source=="UnitSquare")                                     // source is an internally generated dolfin mesh
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitSquareMesh(cells[0], cells[1], diagonal) );

    Side left(0, 0.0);
    Side right(0, 1.0);
    Side bottom(1, 0.0);
    Side top(1, 1.0);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 1, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
    bottom.mark(*edgeids, 3);
    top.mark(*edgeids, 4);
  }
  else if (source=="Rectangle")                                      // source is an internally generated dolfin mesh
  {
    std::vector<double> lowerleft;
    buffer.str(""); buffer << optionpath << "/source/lower_left";
    serr = Spud::get_option(buffer.str(), lowerleft); 
    spud_err(buffer.str(), serr);
    
    std::vector<double> upperright;
    buffer.str(""); buffer << optionpath << "/source/upper_right";
    serr = Spud::get_option(buffer.str(), upperright); 
    spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    std::string diagonal;
    buffer.str(""); buffer << optionpath << "/source/diagonal";
    serr = Spud::get_option(buffer.str(), diagonal); 
    spud_err(buffer.str(), serr);
    
    const dolfin::Point lowerleftpoint(2, lowerleft.data());
    const dolfin::Point upperrightpoint(2, upperright.data());
    mesh.reset(new dolfin::RectangleMesh(lowerleftpoint, upperrightpoint, 
                                    cells[0], cells[1], diagonal));

    Side left(0, lowerleft[0]);
    Side right(0, upperright[0]);
    Side bottom(1, lowerleft[1]);
    Side top(1, upperright[1]);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 1, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
    bottom.mark(*edgeids, 3);
    top.mark(*edgeids, 4);
  }
  else if (source=="UnitCube")                                       // source is an internally generated dolfin mesh
  {
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    mesh.reset( new dolfin::UnitCubeMesh(cells[0], 
                                     cells[1], 
                                     cells[2]) );

    Side left(0, 0.0);
    Side right(0, 1.0);
    Side bottom(2, 0.0);
    Side top(2, 1.0);
    Side back(1, 0.0);
    Side front(1, 1.0);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 2, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
    bottom.mark(*edgeids, 3);
    top.mark(*edgeids, 4);
    back.mark(*edgeids, 5);
    front.mark(*edgeids, 6);
  }
  else if (source=="Box")                                            // source is an internally generated dolfin mesh
  {
    std::vector<double> lowerbackleft;
    buffer.str(""); buffer << optionpath 
                              << "/source/lower_back_left";
    serr = Spud::get_option(buffer.str(), lowerbackleft); 
    spud_err(buffer.str(), serr);
    
    std::vector<double> upperfrontright;
    buffer.str(""); buffer << optionpath 
                            << "/source/upper_front_right";
    serr = Spud::get_option(buffer.str(), upperfrontright); 
    spud_err(buffer.str(), serr);
    
    std::vector<int> cells;
    buffer.str(""); buffer << optionpath << "/source/number_cells";
    serr = Spud::get_option(buffer.str(), cells); 
    spud_err(buffer.str(), serr);
    
    const dolfin::Point lowerbackleftpoint(3, lowerbackleft.data());
    const dolfin::Point upperfrontrightpoint(3, upperfrontright.data());
    mesh.reset( new dolfin::BoxMesh(lowerbackleftpoint, 
                                upperfrontrightpoint, 
                                cells[0], cells[1], cells[2]) );

    Side left(0, lowerbackleft[0]);
    Side right(0, upperfrontright[0]);
    Side bottom(2, lowerbackleft[2]);
    Side top(2, upperfrontright[2]);
    Side back(1, lowerbackleft[1]);
    Side front(1, upperfrontright[1]);

    edgeids.reset(new dolfin::MeshFunction<std::size_t>(mesh, 2, 0));
    left.mark(*edgeids, 1);
    right.mark(*edgeids, 2);
    bottom.mark(*edgeids, 3);
    top.mark(*edgeids, 4);
    back.mark(*edgeids, 5);
    front.mark(*edgeids, 6);
  }
  else                                                               // source is unrecognised
  {
    tf_err("Unknown mesh source.", "Don't understand mesh description.");
  }

  (*mesh).rename(meshname, meshname);
  register_mesh(mesh, meshname, optionpath,
                cellids, edgeids);                                   // put the new mesh in the bucket
}

//*******************************************************************|************************************************************//
// create a new system, fill it and put it into the bucket
//*******************************************************************|************************************************************//
void SpudBucket::fill_systems_(const std::string &optionpath)
{
  SpudSystemBucket_ptr system(new SpudSystemBucket(optionpath,       // create a new system (assumed to be a spudsystem with this 
                                                            this));  // bucket as a parent)

  (*system).fill();                                                  // fill the system

  register_system(system, (*system).name());                         // put the system in the bucket
}

//*******************************************************************|************************************************************//
// loop over the systems defined in the options dictionary and register the base uflsymbols for all coefficient functions
//*******************************************************************|************************************************************//
void SpudBucket::fill_baseuflsymbols_(const std::string &optionpath)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::string uflsymbol, type;                                       // base ufl symbol and coefficient type

  buffer.str("");  buffer << optionpath << "/coefficient";           // the optionpath is for a system so append the coefficient tag
  int ncoeffs = Spud::option_count(buffer.str());                    // and find out how many there are
  for (uint i = 0; i < ncoeffs; i++)                                 // loop over them
  {
    buffer.str(""); buffer << optionpath << "/coefficient[" << i 
                                                    << "]/type/name";
    serr = Spud::get_option(buffer.str(), type);                     // find out the coefficient type
    spud_err(buffer.str(), serr);
    if (type=="Function")                                            // if it's a function then register its derived ufl symbols
    {
      buffer.str(""); buffer << optionpath << "/coefficient[" << i 
                                                   << "]/ufl_symbol";
      serr = Spud::get_option(buffer.str(), uflsymbol); 
      spud_err(buffer.str(), serr);
      register_baseuflsymbol(uflsymbol, uflsymbol);
      register_baseuflsymbol(uflsymbol, uflsymbol+"_i");
      register_baseuflsymbol(uflsymbol, uflsymbol+"_n");
    }
  }
}

//*******************************************************************|************************************************************//
// loop over the detectors defined in the options dictionary and set up the requested detectors
//*******************************************************************|************************************************************//
void SpudBucket::fill_detectors_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  GenericDetectors_ptr det(new GenericDetectors());                  // initialize a pointer to a generic detector

  int npdets = Spud::option_count("/io/detectors/point");            // number of point detectors
  for (uint i=0; i<npdets; i++)                                      // loop over point detectors
  {
    std::stringstream detectorpath;
    detectorpath.str(""); detectorpath << "/io/detectors/point[" 
                                                        << i << "]";
    
    std::string detname;                                             // detector name
    buffer.str(""); buffer << detectorpath.str() << "/name";
    serr = Spud::get_option(buffer.str(), detname);
    spud_err(buffer.str(), serr);

    std::vector<double> point;                                       // point location
    buffer.str(""); buffer << detectorpath.str();
    serr = Spud::get_option(buffer.str(), point);
    spud_err(buffer.str(), serr);
    
    det.reset(new PointDetectors(point, detname));                   // create point detector
    register_detector(det, detname, detectorpath.str());             // register detector
  }
  
  int nadets = Spud::option_count("/io/detectors/array");            // number of array detectors
  for (uint i=0; i<nadets; i++)                                      // loop over array detectors
  {
    std::stringstream detectorpath;
    detectorpath.str(""); detectorpath << "/io/detectors/array[" 
                                                        << i << "]";
    
    std::string detname;                                             // detector array name
    buffer.str(""); buffer << detectorpath.str() << "/name";
    serr = Spud::get_option(buffer.str(), detname);
    spud_err(buffer.str(), serr);

    std::string function;                                            // python function describing the detector positions
    buffer.str(""); buffer << detectorpath.str() << "/python";
    serr = Spud::get_option(buffer.str(), function);
    spud_err(buffer.str(), serr);
    
                                                                     // create python detectors array
    det.reset(new PythonDetectors(dimension(), function, detname));
    register_detector(det, detname, detectorpath.str());             // register detector
  }  
  
}

//*******************************************************************|************************************************************//
// for each solve location check if the user has specified a solution order or not
// if they have read he options, if not, default to using the flags under the solvers
//*******************************************************************|************************************************************//
void SpudBucket::fill_systemssolvers_()
{
  std::stringstream buffer;                                          // optionpath buffer

  std::map<const int, std::string> solve_locations;
  solve_locations[SOLVE_START]       = "at_start";
  solve_locations[SOLVE_TIMELOOP]    = "in_timeloop";
  solve_locations[SOLVE_DIAGNOSTICS] = "with_diagnostics";

  std::map<const int, std::string>::const_iterator sl_it;
  for (sl_it = solve_locations.begin(); sl_it != solve_locations.end(); sl_it++)
  {
    buffer.str(""); buffer << "/solution_order/solve::" << (*sl_it).second;
    if (Spud::have_option(buffer.str()))                               // if the user has specified a solution order
    {                                                                  // then use that
      SpudSystemsSolverBucket_ptr solver( new SpudSystemsSolverBucket(buffer.str(), (*sl_it).first, this) );
      (*solver).fill();

      register_systemssolver(solver, (*sl_it).first);
    }
    else                                                               // otherwise, use a default version that
    {                                                                  // just collects all solvers based on their
      SystemsSolverBucket_ptr solver( new SystemsSolverBucket((*sl_it).first, this) );// individual solve location flags
      (*solver).fill();

      register_systemssolver(solver, (*sl_it).first);
    }
  }

}

//*******************************************************************|************************************************************//
// initialize the data structures for diagnostic output
//*******************************************************************|************************************************************//
void SpudBucket::fill_diagnostics_()
{
  std::stringstream buffer;                                          // optionpath buffer

  statfile_.reset( new StatisticsFile(output_basename()+".stat", 
                           (*(*meshes_begin()).second).mpi_comm(),
                           this) );

  int npdets = Spud::option_count("/io/detectors/point");            // number of point detectors
  int nadets = Spud::option_count("/io/detectors/array");            // number of array detectors
  if ((npdets + nadets) > 0)
  {
    if (Spud::option_count("/system/field/diagnostics/include_in_detectors")==0)
    {
      tf_err("Reqested detectors but selected no fields to include.", 
             "No fields included in detectors.");
    }

    detfile_.reset( new DetectorsFile(output_basename()+".det", 
                           (*(*meshes_begin()).second).mpi_comm(),
                           this) );
  }

  if ((Spud::option_count("/system/functional/include_in_steady_state")+
       Spud::option_count("/system/field/diagnostics/include_in_steady_state"))>0)
  {
    steadyfile_.reset( new SteadyStateFile(output_basename()+".steady",
                           (*(*meshes_begin()).second).mpi_comm(),
                           this) );
  }

  for (SystemBucket_it s_it = systems_begin(); s_it != systems_end(); s_it++)
  {
    (*std::dynamic_pointer_cast< SpudSystemBucket >((*s_it).second)).initialize_diagnostics();
  }

  for (i_SystemsSolverBucket_it ss_it = systemssolvers_begin(); 
                           ss_it != systemssolvers_end(); ss_it++)
  {
    SpudSystemsSolverBucket_ptr sss_ptr = std::dynamic_pointer_cast<SpudSystemsSolverBucket>((*ss_it).second);
    if (sss_ptr)                                                     // only spud system solvers can have diagnostics
    {
      (*sss_ptr).initialize_diagnostics();
    }
  }
}

//*******************************************************************|************************************************************//
// checkpoint the options file
//*******************************************************************|************************************************************//
void SpudBucket::checkpoint_options_(const double_ptr time)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;

  buffer.str(""); buffer << "/io/output_base_name";
  namebuffer.str(""); namebuffer << output_basename() 
                                 << "_checkpoint";
  serr = Spud::set_option(buffer.str(), namebuffer.str());
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << "/timestepping";                         // is this a timestepping simulation?
  if (Spud::have_option(buffer.str()))
  {
    buffer.str(""); buffer << "/timestepping/current_time";          // set the current time
    serr = Spud::set_option(buffer.str(), *time);
    spud_err(buffer.str(), serr);

    buffer.str(""); 
    buffer << "/timestepping/timestep/coefficient::Timestep/type::Constant/rank::Scalar/value::WholeMesh/constant";
    serr = Spud::set_option(buffer.str(), timestep());
    spud_err(buffer.str(), serr);
  }

  if (dolfin::MPI::rank((*(*meshes_begin()).second).mpi_comm())==0)
  {
    namebuffer.str(""); namebuffer << output_basename() 
                                   << "_checkpoint_" 
                                   << checkpoint_count() 
                                   << ".tfml";
    Spud::write_options(namebuffer.str());
  }
  
  buffer.str(""); buffer << "/io/output_base_name";
  serr = Spud::set_option(buffer.str(), output_basename());
  spud_err(buffer.str(), serr);

}


