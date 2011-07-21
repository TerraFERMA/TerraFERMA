
#include "DiagnosticsFile.h"
#include "Bucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>

using namespace buckettools;

DiagnosticsFile::DiagnosticsFile(const std::string name) : StatFile(name)
{
  // Do nothing... all handled by StatFile constructor
}

DiagnosticsFile::~DiagnosticsFile()
{
  // Do nothing... all handled by StatFile destructor
}

void DiagnosticsFile::write_header(const Bucket &bucket, 
                                   const bool &timestepping)
{
  uint column = 1;
  
  file_ << "<header>" << std::endl;
  header_constants_();
  if (timestepping)
  {
    header_timestep_(column);
  }
  header_bucket_(bucket, column);
  file_ << "</header>" << std::endl;
}

void DiagnosticsFile::header_bucket_(const Bucket &bucket, uint &column)
{
  std::stringstream buffer;

  for (SystemBucket_const_it sys_it = bucket.systems_begin(); sys_it != bucket.systems_end(); sys_it++)
  {
    header_functionbucket_((*(*sys_it).second).fields_begin(), (*(*sys_it).second).fields_end(), column);

    header_functionbucket_((*(*sys_it).second).coeffs_begin(), (*(*sys_it).second).coeffs_end(), column);
  }

}

void DiagnosticsFile::header_functionbucket_(FunctionBucket_const_it f_begin, FunctionBucket_const_it f_end, uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    if ((*(*f_it).second).include_in_diagnostics())
    {
      if ((*(*(*f_it).second).function()).value_rank()==0)
      {
        tag_((*(*f_it).second).name(), column, "max", (*(*(*f_it).second).system()).name());
        column++;
        tag_((*(*f_it).second).name(), column, "min", (*(*(*f_it).second).system()).name());
        column++;
      }
      else if ((*(*(*f_it).second).function()).value_rank()==1)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        tag_((*(*f_it).second).name(), column, "max", (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else if ((*(*(*f_it).second).function()).value_rank()==2)
      {
        int components = (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
        tag_((*(*f_it).second).name(), column, "max", (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else
      {
        dolfin::error("In DiagnosticsFile::header_bucket_, unknown function rank.");
      }
    }
  }
}
