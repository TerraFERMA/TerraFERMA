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


#include "SystemsConvergenceFile.h"
#include "Bucket.h"
#include "SystemBucket.h"
#include "SolverBucket.h"
#include "SystemsSolverBucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>
#include "Logger.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SystemsConvergenceFile::SystemsConvergenceFile(const std::string &name, 
                                               const MPI_Comm &comm, 
                                               const Bucket *bucket,
                                               const SystemsSolverBucket *systemssolver) : // specific constructor
                                               DiagnosticsFile(name, comm, bucket), systemssolver_(systemssolver)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SystemsConvergenceFile::~SystemsConvergenceFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::write_data(const double &norm)
{
  initialize_();
  
  data_timestep_();                                                  // write the timestepping information
  data_iteration_();                                                 // write the iteration information
  data_bucket_(norm);                                                // write the bucket data
  
  data_endlineflush_();
  
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::write_header_()
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
void SystemsConvergenceFile::header_iteration_()
{
  
  header_systemssolver_(systemssolver_);
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::header_bucket_()
{
  tag_((*systemssolver_).name(), "res_norm(l2)");

  const std::vector<SolverBucket_ptr>& residualsolvers = 
                                (*systemssolver_).residualsolvers(); // the residual solvers are identified by having as
                                                                     // many unique system names as possible so this will do for us
                                                                     // here too

  std::vector<SolverBucket_ptr>::const_iterator sol_it;
  for (sol_it = residualsolvers.begin(); sol_it != residualsolvers.end(); sol_it++)
  {
    const SystemBucket* p_sys = (**sol_it).system();
    header_system_(p_sys);                                           // write the header for the system itself

    std::pair<const SystemBucket*, std::vector<FunctionBucket_ptr> > fields;
    fields.first = p_sys;
    for (FunctionBucket_const_it f_it = (*p_sys).fields_begin(); 
                                 f_it != (*p_sys).fields_end();      // loop over the given functions
                                 f_it++)
    {
      header_func_((*f_it).second);                                  // write the header for the fields in the system
      fields.second.push_back((*f_it).second);
    }
    fields_.push_back(fields);
  }
}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::header_system_(const SystemBucket* p_sys)
{
  tag_((*p_sys).name(), "max");
  tag_((*p_sys).name(), "min");
  tag_((*p_sys).name(), "res_max");
  tag_((*p_sys).name(), "res_min");
  tag_((*p_sys).name(), "res_norm(l2)");
  tag_((*p_sys).name(), "res_norm(inf)");
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::header_func_(const FunctionBucket_ptr f_ptr)
{
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
  }

}

//*******************************************************************|************************************************************//
// write data to the file for values relating to iterations
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::data_iteration_()
{
  data_systemssolver_(systemssolver_);  
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::data_bucket_(const double &norm)
{
  data_(norm);

  std::vector< std::pair< const SystemBucket*, std::vector<FunctionBucket_ptr> > >::iterator f_it;
  for (f_it = fields_.begin(); 
              f_it != fields_.end(); f_it++)
  {
    data_system_((*f_it).first);

    for (std::vector<FunctionBucket_ptr>::iterator f2_it = (*f_it).second.begin(); 
                                                   f2_it != (*f_it).second.end(); f2_it++)
    {
      data_field_(*f2_it);
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::data_system_(const SystemBucket* p_sys)
{
  std::vector<double> values;

  values.push_back((*(*(*p_sys).iteratedfunction()).vector()).max());
  values.push_back((*(*(*p_sys).iteratedfunction()).vector()).min());
  values.push_back((*(*(*p_sys).residualfunction()).vector()).max());
  values.push_back((*(*(*p_sys).residualfunction()).vector()).min());
  values.push_back((*(*(*p_sys).residualfunction()).vector()).norm("l2"));
  values.push_back((*(*(*p_sys).residualfunction()).vector()).norm("linf"));

  data_(values);
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void SystemsConvergenceFile::data_field_(FunctionBucket_ptr f_ptr)
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
}


