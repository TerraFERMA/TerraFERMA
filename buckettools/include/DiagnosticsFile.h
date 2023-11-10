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


#ifndef __DIAGNOSTICS_FILE_H
#define __DIAGNOSTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include <dolfin.h>

namespace buckettools
{

  class Bucket;
  class SystemsSolverBucket;

  //*****************************************************************|************************************************************//
  // DiagnosticsFile class:
  //
  // This is a base class, which provides simple functionality, for writing various pieces of output to file in a parsable
  // format (using an xml header).  Actual files for specific purposes (diagnostics, detectors etc.) should be defined in
  // classes derived from this base class.
  //*****************************************************************|************************************************************//
  class DiagnosticsFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    DiagnosticsFile(const std::string &name, 
                    const MPI_Comm &comm, 
                    const Bucket *bucket);  // specific constructor
    
    virtual ~DiagnosticsFile();                                      // default destructor
    
    //***************************************************************|***********************************************************//
    // Closing
    //***************************************************************|***********************************************************//

    void close();                                                    // close the file
    
  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::ofstream file_;                                             // file stream

    std::string name_;                                               // file name

    const Bucket *bucket_;                                           // a pointer to the bucket

    const MPI_Comm mpicomm_;                                         // mpi comm

    uint ncolumns_;                                                  // total number of columns

    bool initialized_;

    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void initialize_();

    virtual void write_header_() = 0;

    void header_constants_(const bool &binary=false);                // write the header entries for constant values that do not 
                                                                     // appear later again in the file
    
    void header_open_();                                             // open the header section

    void header_close_();                                            // close the header section

    void header_timestep_();                                         // write the header for a dynamic simulation
                                                                     // (the elapsed time, timestep etc.)
    
    void header_systemssolver_(const SystemsSolverBucket* p_syssol); // write the header for any systems solver iterations

    void constant_tag_(const std::string &name,                      // write a header tag for a constant value
                       const std::string &type, 
                       const std::string &value);
    
    void tag_(const std::string &name,                               // write a header tag for a temporal value
              const std::string &statistic,
              const std::string &system="",
              const uint &components=0);

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void data_endlineflush_();

    void data_timestep_();                                           // write the data for timestepping for a dynamic simulation

    void data_systemssolver_(const SystemsSolverBucket* p_syssol);   // write the data for any systems solver iterations

    void data_(const int &value);                                    // write generic data to the file

    void data_(const double &value);                                 // write generic data to the file

    void data_(const std::vector<double> &values);                   // write generic data to the file
    
  };
  
  typedef std::shared_ptr< DiagnosticsFile > DiagnosticsFile_ptr;  // define a boost shared ptr type for the class

}
#endif
