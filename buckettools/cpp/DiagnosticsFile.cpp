
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

  for (SystemBucket_const_it sys_it = bucket.systems_begin(); sys_it != bucket.systems_end(); sys_it++)
  {
    header_system_(sys_it, column);

    header_field_((*(*sys_it).second).fields_begin(), (*(*sys_it).second).fields_end(), column);

    header_coeff_((*(*sys_it).second).coeffs_begin(), (*(*sys_it).second).coeffs_end(), column);
  }

}

void DiagnosticsFile::header_system_(SystemBucket_const_it sys_it, uint &column)
{
  tag_("SystemFunction", column, "max", (*(*sys_it).second).name());
  column++;
  tag_("SystemFunction", column, "min", (*(*sys_it).second).name());
  column++;
}

void DiagnosticsFile::header_field_(FunctionBucket_const_it f_begin, FunctionBucket_const_it f_end, uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    // loop over the functions (but check that we want to include them in the diagnostic output)
    if ((*(*f_it).second).include_in_diagnostics())
    {
      // if yes, then populate the header with the default stats
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

      header_functional_(f_it, (*(*f_it).second).functionals_begin(), (*(*f_it).second).functionals_end(), column);
    }
  }
}

void DiagnosticsFile::header_coeff_(FunctionBucket_const_it f_begin, FunctionBucket_const_it f_end, uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    // loop over the functions (but check that we want to include them in the diagnostic output)
    if ((*(*f_it).second).include_in_diagnostics())
    {
      header_functional_(f_it, (*(*f_it).second).functionals_begin(), (*(*f_it).second).functionals_end(), column);
    }
  }
}

void DiagnosticsFile::header_functional_(FunctionBucket_const_it f_it, Form_const_it s_begin, Form_const_it s_end, uint &column)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)
  {
    // loop over the functionals attached to this function, adding their names to the header
    tag_((*(*f_it).second).name(), column, (*s_it).first, (*(*(*f_it).second).system()).name());
    column++;
  }
}

void DiagnosticsFile::write_data(Bucket &bucket)
{
  
  data_bucket_(bucket);

  file_ << std::endl << std::flush;
  
}

void DiagnosticsFile::write_data(const uint   &timestep,
                                 const double &elapsedtime, 
                                 const double &dt, 
                                 Bucket       &bucket)
{
  
  data_timestep_(timestep, elapsedtime, dt);
  data_bucket_(bucket);
  
  file_ << std::endl << std::flush;
  
}

void DiagnosticsFile::data_timestep_(const uint   &timestep,
                                     const double &elapsedtime, 
                                     const double &dt)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  file_ << timestep << " ";  
  file_ << elapsedtime << " ";
  file_ << dt << " ";
  
  file_.unsetf(std::ios::scientific);
  
}

void DiagnosticsFile::data_bucket_(Bucket &bucket)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  for (SystemBucket_const_it sys_it = bucket.systems_begin(); sys_it != bucket.systems_end(); sys_it++)
  {
    data_system_(sys_it);

    data_field_((*(*sys_it).second).fields_begin(), (*(*sys_it).second).fields_end());

    data_coeff_((*(*sys_it).second).coeffs_begin(), (*(*sys_it).second).coeffs_end());
  }

  file_.unsetf(std::ios::scientific);
  
}

void DiagnosticsFile::data_system_(SystemBucket_const_it sys_it)
{
  file_ << (*(*(*sys_it).second).function()).vector().max() << " ";
  file_ << (*(*(*sys_it).second).function()).vector().min() << " ";
}

void DiagnosticsFile::data_field_(FunctionBucket_const_it f_begin, FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    // loop over the functions (but check that we want to include them in the diagnostic output)
    if ((*(*f_it).second).include_in_diagnostics())
    {
      // if yes, then populate the header with the default stats
      if ((*(*(*f_it).second).function()).value_rank()==0)
      {
        const Function_ptr func_ptr = boost::dynamic_pointer_cast< dolfin::Function >((*(*f_it).second).function());
        const dolfin::Vector vec = (*func_ptr).vector();
        file_ << vec.max() << " ";
        file_ << vec.min() << " ";
      }
      else if ((*(*(*f_it).second).function()).value_rank()==1)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        for (uint i = 0; i < components; i++)
        {
          file_ << (*boost::dynamic_pointer_cast< dolfin::Function >((*(*f_it).second).function()))[i].vector().max() << " ";
        }
        for (uint i = 0; i < components; i++)
        {
          file_ << (*boost::dynamic_pointer_cast< dolfin::Function >((*(*f_it).second).function()))[i].vector().min() << " ";
        }
      }
      else if ((*(*(*f_it).second).function()).value_rank()==2)
      {
        int dim0 = (*(*(*f_it).second).function()).value_dimension(0);
        int dim1 = (*(*(*f_it).second).function()).value_dimension(1);
        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            file_ << (*boost::dynamic_pointer_cast< dolfin::Function >((*(*f_it).second).function()))[i][j].vector().max() << " ";
          }
        }
        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            file_ << (*boost::dynamic_pointer_cast< dolfin::Function >((*(*f_it).second).function()))[i][j].vector().min() << " ";
          }
        }
      }
      else
      {
        dolfin::error("In DiagnosticsFile::header_bucket_, unknown function rank.");
      }

      data_functional_((*(*f_it).second).functionals_begin(), (*(*f_it).second).functionals_end());
    }
  }
}

void DiagnosticsFile::data_coeff_(FunctionBucket_const_it f_begin, FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    // loop over the functions (but check that we want to include them in the diagnostic output)
    if ((*(*f_it).second).include_in_diagnostics())
    {
      data_functional_((*(*f_it).second).functionals_begin(), (*(*f_it).second).functionals_end());
    }
  }
}

void DiagnosticsFile::data_functional_(Form_const_it s_begin, Form_const_it s_end)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)
  {
    // loop over the functionals attached to this function, assembling them and outputing them to file
    double statistic = dolfin::assemble((*(*s_it).second));
    file_ << statistic << " ";
  }
}

