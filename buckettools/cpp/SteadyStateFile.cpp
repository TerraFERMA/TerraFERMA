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


#include "SteadyStateFile.h"
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
SteadyStateFile::SteadyStateFile(const std::string &name, 
                                 const MPI_Comm &comm, 
                                 const Bucket *bucket) : DiagnosticsFile(name, comm, bucket)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SteadyStateFile::~SteadyStateFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::write_header()
{
  header_open_();
  header_constants_();                                               // write constant tags
  header_timestep_();                                                // write tags for the timesteps
  header_bucket_();                                                  // write tags for the actual bucket variables - fields etc.
  header_close_();
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::write_data()
{
  
  data_timestep_();                                                 // write the timestepping information
  data_bucket_();                                                   // write the bucket data
  
  data_endlineflush_();
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::header_bucket_()
{

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                             sys_it != (*bucket_).systems_end(); 
                             sys_it++)
  {

    for (FunctionBucket_it f_it = (*(*sys_it).second).fields_begin(); 
                           f_it != (*(*sys_it).second).fields_end(); // loop over the given functions
                           f_it++)
    {
      header_func_((*f_it).second);
    }

    for (FunctionBucket_it f_it = (*(*sys_it).second).coeffs_begin(); 
                           f_it != (*(*sys_it).second).coeffs_end(); // loop over the given functions
                           f_it++)
    {
      header_func_((*f_it).second);
    }

  }

}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void SteadyStateFile::header_func_(FunctionBucket_ptr f_ptr)
{

  if ((*f_ptr).include_in_steadystate())                           // check they should be included
  {                                                                // yes, then populate header with steady state change
    if ((*f_ptr).rank()==0)
    {
      tag_((*f_ptr).name(), "change("+((*f_ptr).change_normtype())+")",
                             (*(*f_ptr).system()).name());
    }
    else
    {
      tag_((*f_ptr).name(), "change("+((*f_ptr).change_normtype())+")",
                             (*(*f_ptr).system()).name(), 
                             (*f_ptr).size());
    }
  }
  
  header_functional_(f_ptr);

}

//*******************************************************************|************************************************************//
// write a header for a set of model functionals
//*******************************************************************|************************************************************//
void SteadyStateFile::header_functional_(FunctionBucket_ptr f_ptr)
{
  std::pair<FunctionBucket_ptr, std::vector<Form_const_it> > functionals;
  functionals.first = f_ptr;
  for (Form_const_it s_it = (*f_ptr).functionals_begin(); 
                     s_it != (*f_ptr).functionals_end(); s_it++)     // loop over the functional forms associated with the given
  {                                                                  // function bucket
    if ((*f_ptr).include_functional_in_steadystate((*s_it).first))
    {
      functionals.second.push_back(s_it);
      tag_((*f_ptr).name(), (*s_it).first+"_change",       // write tags for each functional
                                    (*(*f_ptr).system()).name());
    }
  }
  functions_.push_back(functionals);
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::data_bucket_()
{
  
  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin(); 
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
    (*(*sys_it).second).updatechange();
  }

  std::vector< std::pair< FunctionBucket_ptr, std::vector<Form_const_it> > >::iterator f_it;
  for (f_it = functions_.begin(); f_it != functions_.end(); f_it++)
  {
    data_func_((*f_it).first, (*f_it).second);
  }

}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void SteadyStateFile::data_func_(FunctionBucket_ptr f_ptr, 
                                 std::vector<Form_const_it> &functionals)
{
  if ((*f_ptr).changefunction())                                   // check if they should be included in the steady state file
  {                                                                // (shouldn't be true for coefficients)
    const std::size_t lsize = (*f_ptr).size();
    std::vector<double> change(lsize);

    for (uint i = 0; i<lsize; i++)
    {
      change[i] = (*f_ptr).change(i);
    }
    data_(change);
  }

  data_functional_(f_ptr, functionals);
}

//*******************************************************************|************************************************************//
// write data for a set of functional forms
//*******************************************************************|************************************************************//
void SteadyStateFile::data_functional_(FunctionBucket_ptr f_ptr, 
                                       std::vector<Form_const_it> &functionals)
{
  for (std::vector<Form_const_it>::const_iterator s_it = functionals.begin(); 
                                           s_it != functionals.end(); 
                                           s_it++)                   // loop over the given functionals
  {
    const double change = (*f_ptr).functionalchange(*s_it);          // get the value of the functional
    data_(change);
  }
}

