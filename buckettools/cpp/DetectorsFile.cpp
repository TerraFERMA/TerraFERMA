
#include "DetectorsFile.h"
#include "Bucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
DetectorsFile::DetectorsFile(const std::string name) : DiagnosticsFile(name)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
DetectorsFile::~DetectorsFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write header for the model described in the given bucket
//*******************************************************************|************************************************************//
void DetectorsFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  uint column = 1;
  
  file_ << "<header>" << std::endl;
  header_constants_();
  header_timestep_(column);
  header_bucket_(column);
  file_ << "</header>" << std::endl;
}

//*******************************************************************|************************************************************//
// write data for the model described in the attached bucket
//*******************************************************************|************************************************************//
void DetectorsFile::write_data()
{
  
  data_timestep_();
  data_bucket_();
  
  file_ << std::endl << std::flush;
  
}

//*******************************************************************|************************************************************//
// write header for the model described in the given bucket
//*******************************************************************|************************************************************//
void DetectorsFile::header_bucket_(uint &column)
{
  std::stringstream buffer;

  header_detector_((*bucket_).detectors_begin(), 
                   (*bucket_).detectors_end(), 
                   column);

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
    header_func_((*(*sys_it).second).fields_begin(),                 // write the header for the fields in the system
                 (*(*sys_it).second).fields_end(), 
                 (*bucket_).detectors_begin(), 
                 (*bucket_).detectors_end(), column);

    header_func_((*(*sys_it).second).coeffs_begin(),                 // write the header for the coefficients in the system
                 (*(*sys_it).second).coeffs_end(), 
                 (*bucket_).detectors_begin(),
                 (*bucket_).detectors_end(), column);
  }
  
}

//*******************************************************************|************************************************************//
// write header for the detectors
//*******************************************************************|************************************************************//
void DetectorsFile::header_detector_(GenericDetectors_const_it d_begin, 
                                     GenericDetectors_const_it d_end, 
                                     uint &column)
{
  std::stringstream buffer;

  for ( GenericDetectors_const_it d_it = d_begin; d_it != d_end; 
                                                            d_it++ ) // loop over the detectors
  {
    for (uint dim = 0; dim<(*(*d_it).second).dim(); dim++)
    {
      buffer.str("");
      buffer << "position_" << dim;                                  // describe the detector positions
      tag_((*(*d_it).second).name(), column, buffer.str(), 
                                      "", (*(*d_it).second).size());
      column+=(*(*d_it).second).size();
    }
  }
  
}

//*******************************************************************|************************************************************//
// write header for the interaction between the detectors and the given functions
//*******************************************************************|************************************************************//
void DetectorsFile::header_func_(FunctionBucket_const_it f_begin, 
                                 FunctionBucket_const_it f_end, 
                                 GenericDetectors_const_it d_begin,
                                 GenericDetectors_const_it d_end,
                                 uint &column)
{
  std::stringstream buffer;

  for ( FunctionBucket_const_it f_it = f_begin;                      // loop over the functions
                                              f_it != f_end; f_it++)
  {
    
    if ((*(*f_it).second).include_in_detectors())
    {
      for ( GenericDetectors_const_it d_it = d_begin; 
                                    d_it != d_end; d_it++)
      {
        if ((*(*(*f_it).second).function()).value_rank()==0)
        {
          tag_((*(*f_it).second).name(), column, (*(*d_it).second).name(), 
                (*(*(*f_it).second).system()).name(), (*(*d_it).second).size());
          column+=(*(*d_it).second).size();
        }
        else if ((*(*(*f_it).second).function()).value_rank()==1)
        {
          for (uint dim = 0; dim<(*(*(*f_it).second).function()).value_size(); dim++)
          {
            buffer.str("");
            buffer << (*(*f_it).second).name() << "_" << dim;
            tag_(buffer.str(), column, (*(*d_it).second).name(), 
                  (*(*(*f_it).second).system()).name(), (*(*d_it).second).size());
            column+=(*(*d_it).second).size();
          }
        }
        else if ((*(*(*f_it).second).function()).value_rank()==2)
        {
          for (uint dim0 = 0; dim0<(*(*(*f_it).second).function()).value_dimension(0); dim0++)
          {
            for (uint dim1 = 0; dim1<(*(*(*f_it).second).function()).value_dimension(1); dim1++)
            {
              buffer.str("");
              buffer << (*(*f_it).second).name() << "_" << dim0 << "_" << dim1;
              tag_(buffer.str(), column, (*(*d_it).second).name(), 
                   (*(*(*f_it).second).system()).name(), (*(*d_it).second).size());
              column+=(*(*d_it).second).size();
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
  
}

//*******************************************************************|************************************************************//
// write data for the attached bucket
//*******************************************************************|************************************************************//
void DetectorsFile::data_bucket_()
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  data_detector_((*bucket_).detectors_begin(), 
                 (*bucket_).detectors_end());

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
    data_func_((*(*sys_it).second).fields_begin(),                   // write the data for the fields in the system
               (*(*sys_it).second).fields_end(), 
               (*bucket_).detectors_begin(), 
               (*bucket_).detectors_end(),
               (*(*sys_it).second).mesh());

    data_func_((*(*sys_it).second).coeffs_begin(),                   // write the data for the coefficients in the system
               (*(*sys_it).second).coeffs_end(), 
               (*bucket_).detectors_begin(),
               (*bucket_).detectors_end(),
               (*(*sys_it).second).mesh());
  }

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for the detectors
//*******************************************************************|************************************************************//
void DetectorsFile::data_detector_(GenericDetectors_const_it d_begin,
                                   GenericDetectors_const_it d_end)
{
  
  for ( GenericDetectors_const_it d_it = d_begin; 
                                           d_it != d_end; d_it++)
  {
    for (uint dim = 0; dim<(*(*d_it).second).dim(); dim++)
    {
      for (std::vector< Array_double_ptr >::const_iterator pos = 
                                      (*(*d_it).second).begin(); 
                            pos < (*(*d_it).second).end(); pos++)
      {   
        file_ << (**pos)[dim] << " ";
      }
    }
  }
  
}

//*******************************************************************|************************************************************//
// write data for the interaction between the detectors and the functions
//*******************************************************************|************************************************************//
void DetectorsFile::data_func_(FunctionBucket_const_it f_begin, 
                               FunctionBucket_const_it f_end, 
                               GenericDetectors_const_it d_begin,
                               GenericDetectors_const_it d_end,
                               Mesh_ptr mesh)
{
  
  for ( FunctionBucket_const_it f_it = f_begin; f_it != f_end; f_it++)
  {
    
    if ((*(*f_it).second).include_in_detectors())
    {
      for ( GenericDetectors_const_it d_it = d_begin; 
                                             d_it != d_end; d_it++)
      {
        std::vector< Array_double_ptr > values;
        
        GenericFunction_ptr func = (*(*f_it).second).function();

        (*(*d_it).second).eval(values, *func, mesh);
        
        for (uint dim = 0; dim < (*func).value_size(); dim++)
        {
          for(std::vector< Array_double_ptr >::const_iterator val = 
                                                      values.begin();
                                           val < values.end(); val++)
          {
            file_ << (**val)[dim] << " ";
          }
        }
      }
    }
  }
  
}

