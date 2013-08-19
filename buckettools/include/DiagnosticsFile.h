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
  typedef boost::shared_ptr< Bucket > Bucket_ptr;

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

    DiagnosticsFile();                                               // default constructor
    
    DiagnosticsFile(const std::string &name);                        // specific constructor
    
    virtual ~DiagnosticsFile();                                      // default destructor
    
    //***************************************************************|***********************************************************//
    // Closing
    //***************************************************************|***********************************************************//

    const bool is_open() const;                                      // check if file is open

    void close();                                                    // close the file
    
  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // file name

    std::ofstream file_;                                             // file stream

    Bucket_ptr bucket_;                                              // a partial copy of the bucket

    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void header_constants_();                                        // write the header entries for constant values that do not 
                                                                     // appear later again in the file
    
    void header_timestep_(uint &column);                             // write the header for a dynamic simulation (including columns 
                                                                     // for the elapsed time, timestep etc.)
    
    void constant_tag_(const std::string &name,                      // write a header tag for a constant value
                       const std::string &type, 
                       const std::string &value);
    
    void tag_(const std::string &name,                               // write a header tag for a temporal value that does not belong
              const uint &column,                                    // to a system or have components
              const std::string &statistic)
    { tag_(name, column, statistic, "", 0); }

    void tag_(const std::string &name,                               // write a header tag for a temporal value that does not have
              const uint &column,                                    // components
              const std::string &statistic,
              const std::string &system)
    { tag_(name, column, statistic, system, 0); }
    
    void tag_(const std::string &name,                               // write a header tag for a temporal value (most generic)
              const uint &column,
              const std::string &statistic,
              const std::string &system,
              const uint &components);

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void data_timestep_();                                           // write the data for timestepping for a dynamic simulation
    
  };
  
  typedef boost::shared_ptr< DiagnosticsFile > DiagnosticsFile_ptr;  // define a boost shared ptr type for the class

}
#endif
