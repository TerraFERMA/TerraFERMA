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


#ifndef __KSP_CONVERGENCE_FILE_H
#define __KSP_CONVERGENCE_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "DiagnosticsFile.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // predeclarations: a circular dependency between the KSPConvergenceFile class and the Bucket class requires a lot of predeclarations.
  //*****************************************************************|************************************************************//
  class Bucket;
  class SystemBucket;
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;
  class FunctionBucket;
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;
  class SolverBucket;
  typedef boost::shared_ptr< SolverBucket > SolverBucket_ptr;
  typedef std::map< std::string, FunctionBucket_ptr >::const_iterator FunctionBucket_const_it;
  typedef boost::shared_ptr< dolfin::Form > Form_ptr;
  typedef std::map< std::string, Form_ptr >::const_iterator Form_const_it;

  //*****************************************************************|************************************************************//
  // KSPConvergenceFile class:
  //
  // A derived class from the base statfile class intended for the output of diagnostics to file every dump period.
  // Convergence normally include things like norms of residuals.
  //*****************************************************************|************************************************************//
  class KSPConvergenceFile : public DiagnosticsFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    KSPConvergenceFile(const std::string &name, 
                    const std::string &systemname, 
                    const std::string &solvername);                  // specific constructor
 
    ~KSPConvergenceFile();                                              // default destructor
    
    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void write_header(const Bucket &bucket);                         // write header for the bucket

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void write_data(const int &kspit);                               // write data to file for a simulation
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string systemname_;                                         // system name

    std::string solvername_;                                         // solver name

    //***************************************************************|***********************************************************//
    // Header writing functions (continued)
    //***************************************************************|***********************************************************//

    void header_iteration_(uint &column);                            // write the header for iterations 
    
    void header_bucket_(uint &column);                               // write the header for the bucket (non-constant and 
                                                                     // timestepping entries)

    void header_system_(const SystemBucket_ptr sys_ptr,              // write the header for a system
                        uint &column);

    void header_func_(FunctionBucket_const_it f_begin,               // write the header for a set of functions
                      FunctionBucket_const_it f_end, 
                      uint &column);

    //***************************************************************|***********************************************************//
    // Data writing functions (continued)
    //***************************************************************|***********************************************************//

    void data_iteration_(const int &kspit);                          // write the data about the iterations

    void data_bucket_();                                             // write the data for a steady state simulation

    void data_system_(const SystemBucket_ptr sys_ptr);               // write the data for a system

    void data_field_(FunctionBucket_const_it f_begin,                // write the data for a set of fields
                     FunctionBucket_const_it f_end);

  };
  
  typedef boost::shared_ptr< KSPConvergenceFile > KSPConvergenceFile_ptr;    // define a boost shared ptr type for the class

}
#endif
