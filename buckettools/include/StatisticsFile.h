// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __STATISTICS_FILE_H
#define __STATISTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "DiagnosticsFile.h"
#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // predeclarations: a circular dependency between the StatisticsFile class and the Bucket class requires a lot of predeclarations.
  //*****************************************************************|************************************************************//
  class Bucket;
  class SystemBucket;
  typedef std::shared_ptr< SystemBucket > SystemBucket_ptr;
  class FunctionBucket;
  typedef std::shared_ptr< FunctionBucket > FunctionBucket_ptr;
  typedef std::shared_ptr< dolfin::Form > Form_ptr;

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
    
    StatisticsFile(const std::string &name, 
                   const MPI_Comm &comm, 
                   const Bucket *bucket);   // specific constructor
 
    ~StatisticsFile();                                               // default destructor
    
    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void write_header();                         // write header for the bucket

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void write_data();                                               // write data to file for a simulation
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::vector< std::pair< FunctionBucket_ptr, std::vector<Form_const_it> > > functions_;

    //***************************************************************|***********************************************************//
    // Header writing functions (continued)
    //***************************************************************|***********************************************************//

    void header_bucket_();                                           // write the header for the bucket (non-constant and 
                                                                     // timestepping entries)

    void header_func_(const FunctionBucket_ptr f_ptr);               // write the header for a set of functions

    void header_functional_(const FunctionBucket_ptr f_ptr);         // write the header for a set of functionals of a function

    //***************************************************************|***********************************************************//
    // Data writing functions (continued)
    //***************************************************************|***********************************************************//

    void data_bucket_();                                             // write the data for a steady state simulation

    void data_func_(FunctionBucket_ptr f_ptr, 
                    std::vector<Form_const_it> &functionals);        // write the data for a set of functions

    void data_functional_(FunctionBucket_ptr f_ptr,
                          std::vector<Form_const_it> &functionals);  // write the data for a set of functionals

  };
  
  typedef std::shared_ptr< StatisticsFile > StatisticsFile_ptr;    // define a boost shared ptr type for the class

}
#endif
