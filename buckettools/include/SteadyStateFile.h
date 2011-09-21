
#ifndef __STEADYSTATE_FILE_H
#define __STEADYSTATE_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "StatFile.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // SteadyStateFile class:
  //
  // A derived class from the base statfile class intended for the output of diagnostics to file every dump period.
  // Statistics for this file are the change between timesteps.
  //*****************************************************************|************************************************************//
  class SteadyStateFile : public StatFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    SteadyStateFile(const std::string &name);                         // specific constructor
 
    ~SteadyStateFile();                                               // default destructor
    
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
    // Base data
    //***************************************************************|***********************************************************//

    std::string norm_type_;                                          // norm type

    //***************************************************************|***********************************************************//
    // Header writing functions (continued)
    //***************************************************************|***********************************************************//

    void header_bucket_(const Bucket &bucket,                        // write the header for the bucket (non-constant and 
                        uint &column);                               // timestepping entries)

    void header_field_(FunctionBucket_const_it f_begin,              // write the header for a set of fields
                       FunctionBucket_const_it f_end, 
                       uint &column);

    //***************************************************************|***********************************************************//
    // Data writing functions (continued)
    //***************************************************************|***********************************************************//

    void data_bucket_(const Bucket &bucket);                         // write the data for a steady state simulation

    void data_field_(FunctionBucket_const_it f_begin,                // write the data for a set of fields
                     FunctionBucket_const_it f_end);

  };
  
}
#endif
