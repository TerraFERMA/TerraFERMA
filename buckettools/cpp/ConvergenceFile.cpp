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
#include "SystemsSolverBucket.h"
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
                                 const Bucket *bucket,
                                 const std::string &systemname, 
                                 const std::string &solvername,
                                 const std::string &systemssolvername) :
                                      DiagnosticsFile(name, comm, bucket),
                                      systemname_(systemname),
                                      solvername_(solvername),
                                      systemssolvername_(systemssolvername)
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
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_data()
{
  initialize_();
  
  data_timestep_();                                                  // write the timestepping information
  data_iteration_();                                                 // write the iteration information
  data_bucket_();                                                    // write the bucket data
  
  data_endlineflush_();
  
}

//*******************************************************************|************************************************************//
// convenience function to return the associated systems solver (if any)
//*******************************************************************|************************************************************//
SystemsSolverBucket* ConvergenceFile::systemssolver_()
{
  SystemsSolverBucket* p_syssol = NULL;
  if (systemssolvername_.length() > 0)
  {
    p_syssol = (*(*(*bucket_).fetch_system(systemname_)).fetch_solver(solvername_)).fetch_systemssolver(systemssolvername_);
  }
  return p_syssol;
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_header_()
{
  header_open_();
  header_constants_();                                               // write constant tags
  header_timestep_();                                                // write tags for the timesteps
  header_iteration_();                                               // write tags for the iterations
  header_bucket_();                                                  // write tags for the actual bucket variables - fields etc.
  header_close_();
}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::header_iteration_()
{
  SystemsSolverBucket* p_syssol = systemssolver_();
  if (p_syssol)                                                      // any nonlinear systems iterations
  {
    header_systemssolver_(p_syssol);
  }
  tag_("NonlinearIteration", "value");                               // the nonlinear solver iteration
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  header_system_(sys_ptr, sol_ptr);                                  // write the header for the system itself

  for (FunctionBucket_const_it f_it = (*sys_ptr).fields_begin(); 
                               f_it != (*sys_ptr).fields_end();        // loop over the given functions
                               f_it++)
  {
    header_func_((*f_it).second, sol_ptr);
  }

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_system_(const SystemBucket_ptr sys_ptr, 
                                     const SolverBucket_ptr sol_ptr)
{
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
void ConvergenceFile::header_func_(const FunctionBucket_ptr f_ptr, 
                                   const SolverBucket_ptr sol_ptr)
{
  fields_.push_back(f_ptr);

  if ((*f_ptr).rank()==0)             // scalar (no components)
  {
    tag_((*f_ptr).name(), "max", 
                          (*(*f_ptr).system()).name());
    tag_((*f_ptr).name(), "min", 
                          (*(*f_ptr).system()).name());
    tag_((*f_ptr).name(), "res_max", 
                          (*(*f_ptr).system()).name());
    tag_((*f_ptr).name(), "res_min", 
                          (*(*f_ptr).system()).name());
    tag_((*f_ptr).name(), "res_norm(l2)", 
                          (*(*f_ptr).system()).name());
    tag_((*f_ptr).name(), "res_norm(linf)", 
                          (*(*f_ptr).system()).name());
    if ((*sol_ptr).type()=="SNES")
    {
      tag_((*f_ptr).name(), "update_max", 
                            (*(*f_ptr).system()).name());
      tag_((*f_ptr).name(), "update_min", 
                            (*(*f_ptr).system()).name());
    }
  }
  else
  {
    tag_((*f_ptr).name(), "max", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    tag_((*f_ptr).name(), "min", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    tag_((*f_ptr).name(), "res_max", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    tag_((*f_ptr).name(), "res_min", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    tag_((*f_ptr).name(), "res_norm(l2)", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    tag_((*f_ptr).name(), "res_norm(linf)", 
              (*(*f_ptr).system()).name(), 
              (*f_ptr).size());
    if ((*sol_ptr).type()=="SNES")
    {
      tag_((*f_ptr).name(), "update_max", 
                (*(*f_ptr).system()).name(), 
                (*f_ptr).size());
      tag_((*f_ptr).name(), "update_min", 
                (*(*f_ptr).system()).name(), 
                (*f_ptr).size());
    }
  }
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::data_iteration_()
{
  SystemsSolverBucket* p_syssol = systemssolver_();
  if (p_syssol)                                                      // any nonlinear systems iterations
  {
    data_systemssolver_(p_syssol);
  }

  SolverBucket_ptr sol_ptr = (*(*bucket_).fetch_system(systemname_)).fetch_solver(solvername_);
  data_((*sol_ptr).iteration_count());
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::data_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);
  
  data_system_(sys_ptr, sol_ptr);

  for (std::vector<FunctionBucket_ptr>::const_iterator f_it = fields_.begin(); f_it != fields_.end(); f_it++)
  {
    data_field_(*f_it, sol_ptr);
  }

}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void ConvergenceFile::data_system_(const SystemBucket_ptr sys_ptr, 
                                   const SolverBucket_ptr sol_ptr)
{
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

  data_(values);

}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void ConvergenceFile::data_field_(const FunctionBucket_ptr f_ptr, const SolverBucket_ptr sol_ptr) 
{
  const std::size_t lsize = (*f_ptr).size();
  std::vector<double> max(lsize), min(lsize), 
                      l2norm(lsize), linfnorm(lsize);

  for (uint i = 0; i<lsize; i++)
  {
    max[i] = (*f_ptr).max("iterated", i);
    min[i] = (*f_ptr).min("iterated", i);
  }
  data_(max);
  data_(min);

  for (uint i = 0; i<lsize; i++)
  {
    max[i] = (*f_ptr).max("residual", i);
    min[i] = (*f_ptr).min("residual", i);
    l2norm[i]   = (*f_ptr).norm("residual", "l2",   i);
    linfnorm[i] = (*f_ptr).norm("residual", "linf", i);
  }
  data_(max);
  data_(min);
  data_(l2norm);
  data_(linfnorm);

  if ((*sol_ptr).type()=="SNES")
  {
    
    for (uint i = 0; i<lsize; i++)
    {
      max[i] = (*f_ptr).max("snesupdate", i);
      min[i] = (*f_ptr).min("snesupdate", i);
    }
    data_(max);
    data_(min);
  }
}


