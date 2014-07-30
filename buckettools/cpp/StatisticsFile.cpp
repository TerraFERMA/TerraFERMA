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


#include "StatisticsFile.h"
#include "Bucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
StatisticsFile::StatisticsFile(const std::string &name, const MPI_Comm &comm) : DiagnosticsFile(name, comm)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
StatisticsFile::~StatisticsFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  header_open_();
  header_constants_();                                               // write constant tags
  header_timestep_();                                                // write tags for the timesteps
  header_bucket_();                                                  // write tags for the actual bucket variables - fields etc.
  header_close_();
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::write_data()
{
  
  data_timestep_();                                                 // write the timestepping information
  data_bucket_();                                                   // write the bucket data
  
  data_endlineflush_();
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::header_bucket_()
{

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
//    header_system_((*sys_it).second);                               // write the header for the system itself

    header_func_((*(*sys_it).second).fields_begin(),                // write the header for the fields in the system
                          (*(*sys_it).second).fields_end());

    header_func_((*(*sys_it).second).coeffs_begin(),                // write the header for the coefficients in the system
                          (*(*sys_it).second).coeffs_end());
  }

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::header_system_(const SystemBucket_ptr sys_ptr)
{
  tag_((*sys_ptr).name(), "max");
  tag_((*sys_ptr).name(), "min");
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void StatisticsFile::header_func_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given functions
                                                      f_it++)
  {
    if ((*(*f_it).second).include_in_statistics())                   // check they should be included
    {                                                                // yes, then populate header with default stats (min and max)
      if ((*(*f_it).second).rank()==0)
      {
        tag_((*(*f_it).second).name(), "max", 
                               (*(*(*f_it).second).system()).name());
        tag_((*(*f_it).second).name(), "min", 
                               (*(*(*f_it).second).system()).name());
        if ((*(*f_it).second).residualfunction())                      // all fields should get in here
        {
          tag_((*(*f_it).second).name(), "res_max", 
                                 (*(*(*f_it).second).system()).name());
          tag_((*(*f_it).second).name(), "res_min", 
                                 (*(*(*f_it).second).system()).name());
        }
      }
      else
      {
        tag_((*(*f_it).second).name(), "max", 
                               (*(*(*f_it).second).system()).name(), 
                               (*(*f_it).second).size());
        tag_((*(*f_it).second).name(), "min", 
                               (*(*(*f_it).second).system()).name(),
                               (*(*f_it).second).size());
        if ((*(*f_it).second).residualfunction())                      // all fields should get in here
        {
          tag_((*(*f_it).second).name(), "res_max", 
                                 (*(*(*f_it).second).system()).name(),
                                 (*(*f_it).second).size());
          tag_((*(*f_it).second).name(), "res_min", 
                                 (*(*(*f_it).second).system()).name(),
                                 (*(*f_it).second).size());
        }
      }

      header_functional_((*f_it).second, 
                        (*(*f_it).second).functionals_begin(),       // write header for any functionals associated with this field
                        (*(*f_it).second).functionals_end());
    }
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model functionals
//*******************************************************************|************************************************************//
void StatisticsFile::header_functional_(const FunctionBucket_ptr f_ptr, 
                                        Form_const_it f_begin, 
                                        Form_const_it f_end)
{
  for (Form_const_it f_it = f_begin; f_it != f_end; f_it++)          // loop over the functional forms associated with the given
  {                                                                  // function bucket
    tag_((*f_ptr).name(), (*f_it).first,                   // write tags for each functional
                                  (*(*f_ptr).system()).name());
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::data_bucket_()
{
  
  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin(); 
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
//    data_system_((*sys_it).second);

    data_func_((*(*sys_it).second).fields_begin(), 
               (*(*sys_it).second).fields_end());

    data_func_((*(*sys_it).second).coeffs_begin(), 
               (*(*sys_it).second).coeffs_end());

  }

}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void StatisticsFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  double max = (*(*(*sys_ptr).function()).vector()).max();
  data_(max);
  double min = (*(*(*sys_ptr).function()).vector()).min();
  data_(min);
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void StatisticsFile::data_func_(FunctionBucket_const_it f_begin, 
                                FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given functions
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_statistics())                   // check if they should be included in the diagnostics
    {                                                                // yes, start with the default stats... min and max
      const std::size_t lsize = (*(*f_it).second).size();
      std::vector<double> max(lsize), min(lsize);

      for (uint i = 0; i<lsize; i++)
      {
        max[i] = (*(*f_it).second).max("iterated", i);
        min[i] = (*(*f_it).second).min("iterated", i);
      }
      data_(max);
      data_(min);

      if ((*(*f_it).second).residualfunction())                        // all fields should get in here
      {
        for (uint i = 0; i<lsize; i++)
        {
          max[i] = (*(*f_it).second).max("residual", i);
          min[i] = (*(*f_it).second).min("residual", i);
        }
        data_(max);
        data_(min);
      }

      data_functional_((*f_it).second, 
                       (*(*f_it).second).functionals_begin(),        // wttie data for all functionals associated with this field
                            (*(*f_it).second).functionals_end());
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a set of functional forms
//*******************************************************************|************************************************************//
void StatisticsFile::data_functional_(FunctionBucket_ptr f_ptr,
                                       Form_const_it s_begin, 
                                       Form_const_it s_end)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)          // loop over the given functionals
  {
    const double statistic = (*f_ptr).functionalvalue(s_it);         // get the value of the functional
    data_(statistic);
  }
}

