
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
SteadyStateFile::SteadyStateFile(const std::string &name) : StatFile(name)
{
                                                                     // do nothing... all handled by StatFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SteadyStateFile::~SteadyStateFile()
{
                                                                     // do nothing... all handled by StatFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::write_header(const Bucket &bucket)
{
  uint column = 1;                                                   // keep count of how many columns there are
  
  file_ << "<header>" << std::endl;                                  // initialize header xml
  header_constants_();                                               // write constant tags
  header_timestep_(column);                                          // write tags for the timesteps
  header_bucket_(bucket, column);                                    // write tags for the actual bucket variables - fields etc.
  file_ << "</header>" << std::endl << std::flush;                   // finalize header xml
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
// FIXME: this will write data for the bucket in its current state, if this is different to when the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::write_data(const Bucket &bucket)
{
  
  data_timestep_(bucket);                        // write the timestepping information
  data_bucket_(bucket);                                             // write the bucket data
  
  file_ << std::endl << std::flush;                                 // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::header_bucket_(const Bucket &bucket, 
                                     uint &column)
{

  for (SystemBucket_const_it sys_it = bucket.systems_begin();        // loop over the systems
                          sys_it != bucket.systems_end(); sys_it++)
  {
    header_field_((*(*sys_it).second).fields_begin(),                // write the header for the fields in the system
                          (*(*sys_it).second).fields_end(), column);
  }

}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
// FIXME: this will write a header for the bucket in its current state, if this is subsequently modified the header will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::header_field_(FunctionBucket_const_it f_begin, 
                                    FunctionBucket_const_it f_end, 
                                    uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                      f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())                  // check they should be included
    {                                                                // yes, then populate header with stats
      if ((*(*(*f_it).second).function()).value_rank()==0)           // scalar (no components)
      {
        tag_((*(*f_it).second).name(), column++, "change", 
                              (*(*(*f_it).second).system()).name());
      }
      else if ((*(*(*f_it).second).function()).value_rank()==1)      // vector (value_size components)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        tag_((*(*f_it).second).name(), column, "change", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else if ((*(*(*f_it).second).function()).value_rank()==2)      // tensor (value_dimension product components)
      {
        int components = 
          (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
        tag_((*(*f_it).second).name(), column, "change", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
      else                                                           // unknown rank
      {
        dolfin::error("In SteadyStateFile::header_bucket_, unknown function rank.");
      }

    }
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
// FIXME: this will write data for the bucket in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::data_bucket_(const Bucket &bucket)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  for (SystemBucket_const_it sys_it = bucket.systems_begin(); 
                          sys_it != bucket.systems_end(); sys_it++)
  {
    data_field_((*(*sys_it).second).fields_begin(), 
                                (*(*sys_it).second).fields_end());

  }

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for a set of fields
// FIXME: this will write data for the system in its current state, if this has been modified since the header was written the file
// will be invalid
//*******************************************************************|************************************************************//
void SteadyStateFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_steadystate())                  // check if they should be included in the steadystate
    {                                                                // yes, start with the default stats... min and max
      dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).function());
      dolfin::Function oldfunc =                                     // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).oldfunction());
      func.vector() -= oldfunc.vector();                             // a deep copy so this shouldn't affect run values
      if (func.value_rank()==0)                                      // scalars (no components)
      {
        file_ << func.vector().norm(norm_type_) << " ";
      }
      else if (func.value_rank()==1)                                 // vectors (multiple components)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        for (uint i = 0; i < components; i++)
        {
          dolfin::Function funccomp = func[i];                       // take a deep copy of the component of the subfunction
          file_ << funccomp.vector().norm(norm_type_) << " ";
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
            file_ << funccomp.vector().norm(norm_type_) << " ";      // maximum for all components
          }
        }
      }
      else                                                           // unknown rank
      {
        dolfin::error("In SteadyStateFile::data_bucket_, unknown function rank.");
      }

    }
  }
}

