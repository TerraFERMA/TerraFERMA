
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

    header_coeff_((*(*sys_it).second).coeffs_begin(),                // write the header for the coefficients in the system
                          (*(*sys_it).second).coeffs_end(), column);
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

    
    header_functional_((*f_it).second, 
                      (*(*f_it).second).functionals_begin(),         // write header for any functionals associated with this field
                      (*(*f_it).second).functionals_end(), column);
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model coefficients
//*******************************************************************|************************************************************//
void SteadyStateFile::header_coeff_(FunctionBucket_const_it f_begin, 
                                    FunctionBucket_const_it f_end, 
                                    uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                      f_it++)
  {
    header_functional_((*f_it).second, 
                      (*(*f_it).second).functionals_begin(),         // write header for any functionals associated with this coeff
                      (*(*f_it).second).functionals_end(), column);
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model functionals
//*******************************************************************|************************************************************//
void SteadyStateFile::header_functional_(const FunctionBucket_ptr f_ptr, 
                                        Form_const_it f_begin, 
                                        Form_const_it f_end, 
                                        uint &column)
{
  for (Form_const_it f_it = f_begin; f_it != f_end; f_it++)          // loop over the functional forms associated with the given
  {                                                                  // function bucket
    if ((*f_ptr).include_functional_in_steadystate((*f_it).first))
    {
      tag_((*f_ptr).name(), column++, (*f_it).first+"_change",       // write tags for each functional
                                    (*(*f_ptr).system()).name());
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

    data_coeff_((*(*sys_it).second).coeffs_begin(), 
                                (*(*sys_it).second).coeffs_end());
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

    data_functional_((*f_it).second,
                     (*(*f_it).second).functionals_begin(),          // write header for any functionals associated with this field
                     (*(*f_it).second).functionals_end());
  }
}

//*******************************************************************|************************************************************//
// write data for a set of coefficients
//*******************************************************************|************************************************************//
void SteadyStateFile::data_coeff_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    data_functional_((*f_it).second,
                     (*(*f_it).second).functionals_begin(),          // write header for any functionals associated with this coeff
                     (*(*f_it).second).functionals_end());
  }
}

//*******************************************************************|************************************************************//
// write data for a set of functional forms
//*******************************************************************|************************************************************//
void SteadyStateFile::data_functional_(FunctionBucket_ptr f_ptr,
                                       Form_const_it s_begin, 
                                       Form_const_it s_end)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)          // loop over the given functionals
  {
    if ((*f_ptr).include_functional_in_steadystate((*s_it).first))
    {
      double change = (*f_ptr).functionalchange(s_it);               // assemble the functional
      file_ << change << " ";                                     // write to file
    }
  }
}

