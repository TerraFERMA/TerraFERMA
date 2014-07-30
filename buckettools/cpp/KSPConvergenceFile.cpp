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


#include "KSPConvergenceFile.h"
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
KSPConvergenceFile::KSPConvergenceFile(const std::string &name, 
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
KSPConvergenceFile::~KSPConvergenceFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void KSPConvergenceFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  header_open_();
  header_constants_();                                               // write constant tags
  header_timestep_();                                                // write tags for the timesteps
  header_iteration_();                                               // write tags for the iterations
  header_bucket_();                                                  // write tags for the actual bucket variables - fields etc.
  header_close_();
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void KSPConvergenceFile::write_data(const int &kspit)
{
  
  data_timestep_();                                                  // write the timestepping information
  data_iteration_(kspit);                                            // write the iteration information
  data_bucket_();                                                    // write the bucket data
  
  data_endlineflush_();
  
}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to iterations
//*******************************************************************|************************************************************//
void KSPConvergenceFile::header_iteration_()
{
  
  tag_("NonlinearSystemsIteration", "value");                        // the nonlinear systems iteration
  tag_("NonlinearIteration", "value");                               // the nonlinear solver iteration
  tag_("KSPIteration", "value");                                     // the ksp solver iteration
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void KSPConvergenceFile::header_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);

  header_system_(sys_ptr);                                           // write the header for the system itself

  header_func_((*sys_ptr).fields_begin(),                            // write the header for the fields in the system
                        (*sys_ptr).fields_end());

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void KSPConvergenceFile::header_system_(const SystemBucket_ptr sys_ptr)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  tag_((*sys_ptr).name(), "max");
  tag_((*sys_ptr).name(), "min");
  tag_((*sys_ptr).name(), "res_max");
  tag_((*sys_ptr).name(), "res_min");
  tag_((*sys_ptr).name(), "res_norm(l2)");
  tag_((*sys_ptr).name(), "res_norm(inf)");
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void KSPConvergenceFile::header_func_(FunctionBucket_const_it f_begin, 
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
    }
    else
    {
      tag_((*(*f_it).second).name(), "max", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
      tag_((*(*f_it).second).name(), "min", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
      tag_((*(*f_it).second).name(), "res_max", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
      tag_((*(*f_it).second).name(), "res_min", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
      tag_((*(*f_it).second).name(), "res_norm(l2)", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
      tag_((*(*f_it).second).name(), "res_norm(linf)", 
                (*(*(*f_it).second).system()).name(), 
                (*(*f_it).second).size());
    }

  }
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to iterations
//*******************************************************************|************************************************************//
void KSPConvergenceFile::data_iteration_(const int &kspit)
{
  SolverBucket_ptr sol_ptr = (*(*bucket_).fetch_system(systemname_)).fetch_solver(solvername_);
  
  data_((*bucket_).iteration_count());  
  data_((*sol_ptr).iteration_count());
  data_(kspit);
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void KSPConvergenceFile::data_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);
  
  data_system_(sys_ptr);
  data_field_((*sys_ptr).fields_begin(), (*sys_ptr).fields_end());

}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void KSPConvergenceFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  std::vector<double> values;

  values.push_back((*(*(*sys_ptr).iteratedfunction()).vector()).max());
  values.push_back((*(*(*sys_ptr).iteratedfunction()).vector()).min());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).max());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).min());
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).norm("l2"));
  values.push_back((*(*(*sys_ptr).residualfunction()).vector()).norm("linf"));

  data_(values);
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void KSPConvergenceFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  SolverBucket_ptr sol_ptr = (*(*(*f_begin).second).system()).fetch_solver(solvername_);

  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    
    const std::size_t lsize = (*(*f_it).second).size();
    std::vector<double> max(lsize), min(lsize), 
                        l2norm(lsize), linfnorm(lsize);

    for (uint i = 0; i<lsize; i++)
    {
      max[i] = (*(*f_it).second).max("iterated", i);
      min[i] = (*(*f_it).second).min("iterated", i);
    }
    data_(max);
    data_(min);

    for (uint i = 0; i<lsize; i++)
    {
      max[i] = (*(*f_it).second).max("residual", i);
      min[i] = (*(*f_it).second).min("residual", i);
      l2norm[i]   = (*(*f_it).second).norm("residual", "l2",   i);
      linfnorm[i] = (*(*f_it).second).norm("residual", "linf", i);
    }
    data_(max);
    data_(min);
    data_(l2norm);
    data_(linfnorm);
  }
}


