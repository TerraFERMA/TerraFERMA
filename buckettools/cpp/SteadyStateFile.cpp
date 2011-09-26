
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
SteadyStateFile::SteadyStateFile(const std::string &name) : DiagnosticsFile(name)
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
void SteadyStateFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  uint column = 1;                                                   // keep count of how many columns there are
  
  file_ << "<header>" << std::endl;                                  // initialize header xml
  header_constants_();                                               // write constant tags
  header_timestep_(column);                                          // write tags for the timesteps
  header_bucket_(column);                                            // write tags for the actual bucket variables - fields etc.
  file_ << "</header>" << std::endl << std::flush;                   // finalize header xml
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::write_data()
{
  
  data_timestep_();                                                 // write the timestepping information
  data_bucket_();                                                   // write the bucket data
  
  file_ << std::endl << std::flush;                                 // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::header_bucket_(uint &column)
{

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
    header_field_((*(*sys_it).second).fields_begin(),                // write the header for the fields in the system
                          (*(*sys_it).second).fields_end(), column);
  }

}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void SteadyStateFile::header_field_(FunctionBucket_const_it f_begin, 
                                    FunctionBucket_const_it f_end, 
                                    uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                      f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())                  // check they should be included
    {                                                                // yes, then populate header with steady state change
      tag_((*(*f_it).second).name(), column++, "change("+((*(*f_it).second).change_normtype())+")", 
                            (*(*(*f_it).second).system()).name());
    }
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void SteadyStateFile::data_bucket_()
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin(); 
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
    (*(*sys_it).second).updatechange();
    data_field_((*(*sys_it).second).fields_begin(), 
                                (*(*sys_it).second).fields_end());
  }

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void SteadyStateFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())                  // check if they should be included in the steady state file
    {
      file_ << (*(*f_it).second).change() << " ";
    }
  }
}

