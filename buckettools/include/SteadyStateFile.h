
#ifndef __STEADYSTATE_FILE_H
#define __STEADYSTATE_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "DiagnosticsFile.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // predeclarations: a circular dependency between the SteadyStateFile class and the Bucket class requires a lot of predeclarations.
  //*****************************************************************|************************************************************//
  class Bucket;
  class FunctionBucket;
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;
  typedef boost::shared_ptr< dolfin::Form > Form_ptr;
  typedef std::map< std::string, FunctionBucket_ptr >::const_iterator FunctionBucket_const_it;
  typedef std::map< std::string, Form_ptr >::const_iterator Form_const_it;

  //*****************************************************************|************************************************************//
  // SteadyStateFile class:
  //
  // A derived class from the base statfile class intended for the output of diagnostics to file every dump period.
  // Statistics normally include things like function mins and maxes as well as functional output.
  //*****************************************************************|************************************************************//
  class SteadyStateFile : public DiagnosticsFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    SteadyStateFile(const std::string &name);                        // specific constructor
 
    ~SteadyStateFile();                                              // default destructor
    
    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void write_header(const Bucket &bucket);                         // write header for the bucket

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void write_data();                                               // write data to file for a simulation
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class

    //***************************************************************|***********************************************************//
    // Header writing functions (continued)
    //***************************************************************|***********************************************************//

    void header_bucket_(uint &column);                               // write the header for the bucket (non-constant and 
                                                                     // timestepping entries)

    void header_field_(FunctionBucket_const_it f_begin,              // write the header for a set of fields
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_coeff_(FunctionBucket_const_it f_begin,              // write the header for a set of coefficients
                       FunctionBucket_const_it f_end, 
                       uint &column);

    void header_functional_(const FunctionBucket_ptr f_ptr,          // write the header for a set of functionals
                            Form_const_it f_begin, 
                            Form_const_it f_end, 
                            uint &column);

    //***************************************************************|***********************************************************//
    // Data writing functions (continued)
    //***************************************************************|***********************************************************//

    void data_bucket_();                                             // write the data for a steady state simulation

    void data_field_(FunctionBucket_const_it f_begin,                // write the data for a set of fields
                     FunctionBucket_const_it f_end);

    void data_coeff_(FunctionBucket_const_it f_begin,                // write the data for a set of coefficients
                     FunctionBucket_const_it f_end);

    void data_functional_(FunctionBucket_ptr f_ptr,                  // write the data for a set of functionals
                          Form_const_it f_begin, 
                          Form_const_it f_end);
  };
  
  typedef boost::shared_ptr< SteadyStateFile > SteadyStateFile_ptr;  // define a boost shared ptr type for the class

}
#endif
