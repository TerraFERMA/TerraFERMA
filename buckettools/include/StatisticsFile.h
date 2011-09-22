
#ifndef __STATISTICS_FILE_H
#define __STATISTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "DiagnosticsFile.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // StatisticsFile class:
  //
  // A derived class from the base statfile class intended for the output of diagnostics to file every dump period.
  // Statistics normally include things like function mins and maxes as well as functional output.
  //*****************************************************************|************************************************************//
  class StatisticsFile : public DiagnosticsFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    StatisticsFile(const std::string &name);                         // specific constructor
 
    ~StatisticsFile();                                               // default destructor
    
    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void write_header(const Bucket &bucket);                         // write header for the bucket

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void write_data(const Bucket &bucket);                           // write data to file for a simulation
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class

    //***************************************************************|***********************************************************//
    // Header writing functions (continued)
    //***************************************************************|***********************************************************//

    void header_bucket_(const Bucket &bucket,                        // write the header for the bucket (non-constant and 
                        uint &column);                               // timestepping entries)

    void header_system_(const SystemBucket_ptr sys_ptr,              // write the header for a system
                        uint &column);

    void header_field_(FunctionBucket_const_it f_begin,              // write the header for a set of fields
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_coeff_(FunctionBucket_const_it f_begin,              // write the header for a set of coefficients
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_functional_(const FunctionBucket_ptr f_ptr,          // write the header for a set of functionals of a function
                            Form_const_it f_begin,
                            Form_const_it f_end, 
                            uint &column);

    //***************************************************************|***********************************************************//
    // Data writing functions (continued)
    //***************************************************************|***********************************************************//

    void data_bucket_(const Bucket &bucket);                         // write the data for a steady state simulation

    void data_system_(const SystemBucket_ptr sys_ptr);               // write the data for a system

    void data_field_(FunctionBucket_const_it f_begin,                // write the data for a set of fields
                     FunctionBucket_const_it f_end);

    void data_coeff_(FunctionBucket_const_it f_begin,                // write the data for a set of coefficients
                     FunctionBucket_const_it f_end);

    void data_functional_(Form_const_it s_begin,                     // write the data for a set of functionals
                          Form_const_it s_end);

  };
  
  typedef boost::shared_ptr< StatisticsFile > StatisticsFile_ptr;    // define a boost shared ptr type for the class

}
#endif
