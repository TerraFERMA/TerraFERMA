
#include "DiagnosticsFile.h"
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
DiagnosticsFile::DiagnosticsFile(const std::string name) : StatFile(name)
{
                                                                     // do nothing... all handled by StatFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
DiagnosticsFile::~DiagnosticsFile()
{
                                                                     // do nothing... all handled by StatFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket, indicating if this is a dynamic simulation or not
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::write_header(const Bucket &bucket, 
                                   const bool &timestepping)
{
  uint column = 1;                                                   // keep count of how many columns there are
  
  file_ << "<header>" << std::endl;                                  // initialize header xml
  header_constants_();                                               // write constant tags
  if (timestepping)                                                  // if this is a timestepping simulation
  {
    header_timestep_(column);                                        // write tags for the timesteps
  }
  header_bucket_(bucket, column);                                    // write tags for the actual bucket variables - fields etc.
  file_ << "</header>" << std::endl;                                 // finalize header xml
}

//*******************************************************************|************************************************************//
// write data for the (presumed steady state) model described in the given bucket
// FIXME: this will write data for the bucket in its current state, if this is different to when the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::write_data(Bucket &bucket)
{
  
  data_bucket_(bucket);                                              // write the bucket data

  file_ << std::endl << std::flush;                                  // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write data for the (presumed dynamic) model described in the given bucket
// FIXME: this will write data for the bucket in its current state, if this is different to when the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::write_data(const uint   &timestep,
                                 const double &elapsedtime, 
                                 const double &dt, 
                                 Bucket       &bucket)
{
  
  data_timestep_(timestep, elapsedtime, dt);                        // write the timestepping information
  data_bucket_(bucket);                                             // write the bucket data
  
  file_ << std::endl << std::flush;                                 // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_bucket_(const Bucket &bucket, 
                                     uint &column)
{

  for (SystemBucket_const_it sys_it = bucket.systems_begin();        // loop over the systems
                          sys_it != bucket.systems_end(); sys_it++)
  {
    header_system_((*sys_it).second, column);                        // write the header for the system itself

    header_field_((*(*sys_it).second).fields_begin(),                // write the header for the fields in the system
                          (*(*sys_it).second).fields_end(), column);

    header_coeff_((*(*sys_it).second).coeffs_begin(),                // write the header for the coefficients in the system
                          (*(*sys_it).second).coeffs_end(), column);
  }

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_system_(const SystemBucket_ptr sys_ptr, 
                                     uint &column)
{
  tag_((*sys_ptr).name(), column++, "max");
  tag_((*sys_ptr).name(), column++, "min");
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_field_(FunctionBucket_const_it f_begin, 
                                    FunctionBucket_const_it f_end, 
                                    uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                      f_it++)
  {
    if ((*(*f_it).second).include_in_diagnostics())                  // check they should be included
    {                                                                // yes, then populate header with default stats (min and max)
      if ((*(*(*f_it).second).function()).value_rank()==0)           // scalar (no components)
      {
        tag_((*(*f_it).second).name(), column++, "max", 
                              (*(*(*f_it).second).system()).name());
        tag_((*(*f_it).second).name(), column++, "min", 
                              (*(*(*f_it).second).system()).name());
      }
      else if ((*(*(*f_it).second).function()).value_rank()==1)      // vector (value_size components)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        tag_((*(*f_it).second).name(), column, "max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else if ((*(*(*f_it).second).function()).value_rank()==2)      // tensor (value_dimension product components)
      {
        int components = 
          (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
        tag_((*(*f_it).second).name(), column, "max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else                                                           // unknown rank
      {
        dolfin::error("In DiagnosticsFile::header_bucket_, unknown function rank.");
      }

      header_functional_((*f_it).second, 
                        (*(*f_it).second).functionals_begin(),       // write header for any functionals associated with this field
                        (*(*f_it).second).functionals_end(), column);
    }
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model coefficients
// these don't get automatic min and max as they're user prescribed (also only coefficient functions have a vector on which they
// could be easily evaluated)
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_coeff_(FunctionBucket_const_it f_begin, 
                                    FunctionBucket_const_it f_end, 
                                    uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given coefficients
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_diagnostics())                  // check if they are to be included or not
    {
      header_functional_((*f_it).second, 
                        (*(*f_it).second).functionals_begin(),       // write a header for all the functionals associated with this
                        (*(*f_it).second).functionals_end(), column);// coefficient
    }
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model functionals
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_functional_(const FunctionBucket_ptr f_ptr, 
                                        Form_const_it f_begin, 
                                        Form_const_it f_end, 
                                        uint &column)
{
  for (Form_const_it f_it = f_begin; f_it != f_end; f_it++)          // loop over the functional forms associated with the given
  {                                                                  // function bucket
    tag_((*f_ptr).name(), column++, (*f_it).first,                   // write tags for each functional
                                  (*(*f_ptr).system()).name());
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
// FIXME: this will write data for the bucket in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_bucket_(Bucket &bucket)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  for (SystemBucket_const_it sys_it = bucket.systems_begin(); 
                          sys_it != bucket.systems_end(); sys_it++)
  {
    data_system_((*sys_it).second);

    data_field_((*(*sys_it).second).fields_begin(), 
                                (*(*sys_it).second).fields_end());

    data_coeff_((*(*sys_it).second).coeffs_begin(), 
                                (*(*sys_it).second).coeffs_end());
  }

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for a system
// FIXME: this will write data for the system in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  file_ << (*(*sys_ptr).function()).vector().max() << " ";
  file_ << (*(*sys_ptr).function()).vector().min() << " ";
}

//*******************************************************************|************************************************************//
// write data for a set of fields
// FIXME: this will write data for the system in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_diagnostics())                  // check if they should be included in the diagnostics
    {                                                                // yes, start with the default stats... min and max
      dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).function());
      if (func.value_rank()==0)                                      // scalars (no components)
      {
        file_ << func.vector().max() << " ";
        file_ << func.vector().min() << " ";
      }
      else if (func.value_rank()==1)                                 // vectors (multiple components)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        for (uint i = 0; i < components; i++)
        {
          dolfin::Function funccomp = func[i];                       // take a deep copy of the component of the subfunction
          file_ << funccomp.vector().max() << " ";                   // maximum for all components
        }
        for (uint i = 0; i < components; i++)
        {
          dolfin::Function funccomp = func[i];                       // take a deep copy of the component of the subfunction
          file_ << funccomp.vector().min() << " ";                   // minimum for all components
        }
      }
      else if (func.value_rank()==2)                                 // tensor (multiple components)
      {
        int dim0 = (*(*(*f_it).second).function()).value_dimension(0);
        int dim1 = (*(*(*f_it).second).function()).value_dimension(1);
        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function funccomp = func[i][j];                  // take a deep copy of the ijth component of the subfunction
            file_ << funccomp.vector().max() << " ";                 // maximum for all components
          }
        }
        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function funccomp = func[i][j];                  // take a deep copy of the ijth component of the subfunction
            file_ << funccomp.vector().min() << " ";                 // minimum for all components
          }
        }
      }
      else                                                           // unknown rank
      {
        dolfin::error("In DiagnosticsFile::data_bucket_, unknown function rank.");
      }

      data_functional_((*(*f_it).second).functionals_begin(),        // wttie data for all functionals associated with this field
                            (*(*f_it).second).functionals_end());
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a set of coefficients
// FIXME: this will write data for the system in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_coeff_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given coefficients
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_diagnostics())                  // check if this coefficient is to be included in diagnostics
    {
      data_functional_((*(*f_it).second).functionals_begin(),        // write data for all functionals associated with this
                              (*(*f_it).second).functionals_end());  // coefficient
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a set of functional forms
// FIXME: this will write data for the system in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_functional_(Form_const_it s_begin, 
                                       Form_const_it s_end)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)          // loop over the given functionals
  {
    double statistic = dolfin::assemble((*(*s_it).second));          // assemble the functional
    file_ << statistic << " ";                                       // write to file
  }
}

