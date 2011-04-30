
#include "DetectorsFile.h"
#include "Detectors.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>

using namespace buckettools;

DetectorsFile::DetectorsFile(const std::string name) : StatFile(name)
{
  // Do nothing... all handled by StatFile constructor
}

void DetectorsFile::write_header(Bucket &bucket, 
                                 const bool timestepping)
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

void DetectorsFile::header_timestep_(uint &column)
{
  
  tag_("timestep", column++, "value");
  tag_("ElapsedTime", column++, "value");
  tag_("dt", column++, "value");
  
}

void DetectorsFile::header_bucket_(Bucket &bucket, uint &column)
{
  std::stringstream buffer;

  // the detector positions
  for ( std::vector< Detectors_ptr >::iterator det = bucket.detectors_begin(); 
                                           det < bucket.detectors_end(); det++)
  {
    for (uint dim = 0; dim<(**det).dim(); dim++)
    {
      buffer.str("");
      buffer << "position_" << dim;
      tag_((**det).name(), column++, buffer.str(), (**det).size());
    }
  }
  
  // the functions interacting with the detectors
  for ( std::vector< GenericFunction_ptr >::iterator func = bucket.functions_begin(); 
                                                func < bucket.functions_end(); func++)
  {
    
    for ( std::vector< Detectors_ptr >::iterator det = bucket.detectors_begin(); 
                                            det < bucket.detectors_end(); det++)
    {
      if ((**func).value_rank()==0)
      {
        tag_((**func).name(), column++, (**det).name(), (**det).size());
      }
      else if ((**func).value_rank()==1)
      {
        for (uint dim = 0; dim<(**func).value_size(); dim++)
        {
          buffer.str("");
          buffer << (**func).name() << "_" << dim;
          tag_(buffer.str(), column++, (**det).name(), (**det).size());
        }
      }
      else if ((**func).value_rank()==2)
      {
        for (uint dim0 = 0; dim0<(**func).value_dimension(0); dim0++)
        {
          for (uint dim1 = 0; dim1<(**func).value_dimension(1); dim1++)
          {
            buffer.str("");
            buffer << (**func).name() << "_" << dim0 << "_" << dim1;
            tag_(buffer.str(), column++, (**det).name(), (**det).size());
          }
        }
      }
      else
      {
        dolfin::error("In DetectorsFile::header_detectors_, unknown function rank.");
      }
    }
  }
  
}

void DetectorsFile::write_data(Bucket &bucket)
{
  
  data_bucket_(bucket);

  file_ << std::endl << std::flush;
  
}

void DetectorsFile::write_data(const uint   timestep,
                               const double elapsedtime, 
                               const double dt, 
                               Bucket       &bucket)
{
  
  data_timestep_(timestep, elapsedtime, dt);
  data_bucket_(bucket);
  
  file_ << std::endl << std::flush;
  
}

void DetectorsFile::data_timestep_(const uint   timestep,
                                         const double elapsedtime, 
                                         const double dt)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  file_ << timestep << " ";  
  file_ << elapsedtime << " ";
  file_ << dt << " ";
  
  file_.unsetf(std::ios::scientific);
  
}

void DetectorsFile::data_bucket_(Bucket &bucket)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  // the detector positions
  for ( std::vector< Detectors_ptr >::iterator det = bucket.detectors_begin(); 
                                           det < bucket.detectors_end(); det++)
  {
    for (uint dim = 0; dim<(**det).dim(); dim++)
    {
      for (std::vector< Array_double_ptr >::iterator pos = (**det).begin();
                                                pos < (**det).end(); pos++)
      {   
        file_ << (**pos)[dim] << " ";
      }
    }
  }
  
  // the functions interacting with the detectors
  for ( std::vector< GenericFunction_ptr >::iterator func = bucket.functions_begin(); 
                                                func < bucket.functions_end(); func++)
  {
    
    for ( std::vector< Detectors_ptr >::iterator det = bucket.detectors_begin(); 
                                            det < bucket.detectors_end(); det++)
    {
      std::vector< Array_double_ptr > values;
      
      (**det).eval(values, **func);
      
      for (uint dim = 0; dim<(**func).value_size(); dim++)
      {
        for(std::vector< Array_double_ptr >::iterator val = values.begin();
                                              val < values.end(); val++)
        {
          file_ << (**val)[dim] << " ";
        }
      }      
    }
  }
    
  file_.unsetf(std::ios::scientific);
  
}

