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


#include "ConvergenceFile.h"
#include "Bucket.h"
#include "SystemBucket.h"
#include "SolverBucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
ConvergenceFile::ConvergenceFile(const std::string &name, 
                                 const MPI_Comm &comm,
                                 const std::string &systemname, 
                                 const std::string &solvername) :
                                      DiagnosticsFile(name, comm),
                                      systemname_(systemname),
                                      solvername_(solvername)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
ConvergenceFile::~ConvergenceFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "<header>" << std::endl;                                // initialize header xml
  }
  header_constants_();                                               // write constant tags
  header_timestep_();                                                // write tags for the timesteps
  header_iteration_();                                               // write tags for the iterations
  header_bucket_();                                                  // write tags for the actual bucket variables - fields etc.
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "</header>" << std::endl << std::flush;                 // finalize header xml
  }
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_data()
{
  
  data_timestep_();                                                  // write the timestepping information
  data_iteration_();                                                 // write the iteration information
  data_bucket_();                                                    // write the bucket data
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << std::endl << std::flush;                                // flush the buffer
  }
  
}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::header_iteration_()
{
  
  tag_("NonlinearSystemsIteration", "value");                        // the nonlinear systems iteration
  tag_("NonlinearIteration", "value");                               // the nonlinear solver iteration
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);

  header_system_(sys_ptr);                                           // write the header for the system itself

  header_func_((*sys_ptr).fields_begin(),                            // write the header for the fields in the system
                        (*sys_ptr).fields_end());

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_system_(const SystemBucket_ptr sys_ptr)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  tag_((*sys_ptr).name(), "max");
  tag_((*sys_ptr).name(), "min");
  tag_((*sys_ptr).name(), "res_max");
  tag_((*sys_ptr).name(), "res_min");
  tag_((*sys_ptr).name(), "res_norm(l2)");
  tag_((*sys_ptr).name(), "res_norm(inf)");
  if ((*sol_ptr).type()=="SNES")
  {
    tag_((*sys_ptr).name(), "update_max");
    tag_((*sys_ptr).name(), "update_min");
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void ConvergenceFile::header_func_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  SolverBucket_ptr sol_ptr = (*(*(*f_begin).second).system()).fetch_solver(solvername_);

  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given functions
                                                      f_it++)
  {
    if ((*(*(*f_it).second).function()).value_rank()==0)             // scalar (no components)
    {
      tag_((*(*f_it).second).name(), "max", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), "min", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), "res_max", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), "res_min", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), "res_norm(l2)", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), "res_norm(linf)", 
                            (*(*(*f_it).second).system()).name());
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), "update_max", 
                              (*(*(*f_it).second).system()).name());
        tag_((*(*f_it).second).name(), "update_min", 
                              (*(*(*f_it).second).system()).name());
      }
    }
    else if ((*(*(*f_it).second).function()).value_rank()==1)        // vector (value_size components)
    {
      int components = (*(*(*f_it).second).function()).value_size();
      tag_((*(*f_it).second).name(), "max", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "min", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_max", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_min", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_norm(l2)", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_norm(linf)", 
                (*(*(*f_it).second).system()).name(), components);
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), "update_max", 
                  (*(*(*f_it).second).system()).name(), components);
        tag_((*(*f_it).second).name(), "update_min", 
                  (*(*(*f_it).second).system()).name(), components);
      }
    }
    else if ((*(*(*f_it).second).function()).value_rank()==2)        // tensor (value_dimension product components)
    {
      int components = 
        (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
      tag_((*(*f_it).second).name(), "max", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "min", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_max", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_min", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_norm(l2)", 
                (*(*(*f_it).second).system()).name(), components);
      tag_((*(*f_it).second).name(), "res_norm(linf)", 
                (*(*(*f_it).second).system()).name(), components);
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), "update_max", 
                  (*(*(*f_it).second).system()).name(), components);
        tag_((*(*f_it).second).name(), "update_min", 
                  (*(*(*f_it).second).system()).name(), components);
      }
    }
    else                                                             // unknown rank
    {
      dolfin::error("In ConvergenceFile::header_bucket_, unknown function rank.");
    }

  }
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::data_iteration_()
{
  SolverBucket_ptr sol_ptr = (*(*bucket_).fetch_system(systemname_)).fetch_solver(solvername_);
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << (*bucket_).iteration_count() << " ";  
    file_ << (*sol_ptr).iteration_count() << " ";
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::data_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.setf(std::ios::scientific);
    file_.precision(10);
  }
  
  data_system_(sys_ptr);

  data_field_((*sys_ptr).fields_begin(), (*sys_ptr).fields_end());

  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.unsetf(std::ios::scientific);
  }
  
}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void ConvergenceFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  std::vector<double> values;

  values.push_back((*(*(*sys_ptr).iteratedfunction()).vector()).max());
  values.push_back((*(*(*sys_ptr).iteratedfunction()).vector()).min());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).max());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).min());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).norm("l2"));
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).norm("linf"));
  if ((*sol_ptr).type()=="SNES")
  {
    values.push_back((*(*(*sys_ptr).snesupdatefunction()).vector()).max());
    values.push_back((*(*(*sys_ptr).snesupdatefunction()).vector()).min());
  }

  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    for (std::vector<double>::const_iterator v = values.begin(); 
                                                  v < values.end(); v++)
    {
      file_ << *v << " ";
    }
  }

}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void ConvergenceFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  SolverBucket_ptr sol_ptr = (*(*(*f_begin).second).system()).fetch_solver(solvername_);

  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    dolfin::Function func =                                          // take a deep copy of the subfunction so the vector is accessible
      *std::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).iteratedfunction());

    data_function_(func);

    dolfin::Function resfunc =                                       // take a deep copy of the subfunction so the vector is accessible
      *std::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).residualfunction());

    data_function_(resfunc, true);

    if ((*sol_ptr).type()=="SNES")
    {
      dolfin::Function upfunc =                                      // take a deep copy of the subfunction so the vector is accessible
        *std::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).snesupdatefunction());

      data_function_(upfunc);
    }
  }
}


