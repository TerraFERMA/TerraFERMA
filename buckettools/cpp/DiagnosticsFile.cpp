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


#include "Bucket.h"
#include "DiagnosticsFile.h"
#include "builddefs.h"
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
DiagnosticsFile::DiagnosticsFile()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
DiagnosticsFile::DiagnosticsFile(const std::string &name) : name_(name)
{
  file_.open((char*)name.c_str());                                   // open the file_ member
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
DiagnosticsFile::~DiagnosticsFile()
{
  close();                                                           // close the file_ member
}

//*******************************************************************|************************************************************//
// write lines of the xml header for constants that do not vary throughout a simulation
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_constants_()
{
  std::stringstream buffer;
  const int maxlencbuffer=1024;
  char cbuffer[maxlencbuffer];
  int cerr;
  
  constant_tag_("GitHash", "string", __GIT_SHA__);                   // the git sha
  
  buffer.str("");
  buffer << __DATE__;
  buffer << " ";
  buffer << __TIME__;
  constant_tag_("CompileTime", "string", buffer.str());              // the comilation time
  
  buffer.str("");
  buffer << ctime(((*bucket_).start_walltime()));
  constant_tag_("StartTime", "string", 
                buffer.str().substr( 0, buffer.str().length() - 1)); // the simulation start time
  
  buffer.str("");
  cerr = gethostname(cbuffer, maxlencbuffer);
  if(cerr==0)
  {
    buffer << cbuffer;
  }
  else
  {
    buffer << "not_defined";
  }
  constant_tag_("HostName", "string", buffer.str());                 // the hostname

}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to timestepping
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_timestep_(uint &column)
{
  
  tag_("timestep", column++, "value");                               // the timestep count (i.e. number of timesteps taken)
  tag_("ElapsedTime", column++, "value");                            // the current time
  tag_("ElapsedWallTime", column++, "value");                        // the elapsed wall time
  tag_("dt", column++, "value");                                     // the actual timestep
  
}

//*******************************************************************|************************************************************//
// write an xml tag for a constant that does not vary throughout a simulation
//*******************************************************************|************************************************************//
void DiagnosticsFile::constant_tag_(const std::string &name, 
                             const std::string &type, 
                             const std::string &value)
{
  
  file_ << "<constant name=\"" << name
        << "\" type=\"" << type
        << "\" value=\"" << value << "\" />" 
        << std::endl << std::flush;

}

//*******************************************************************|************************************************************//
// write an xml tag for a variable in a simulation 
//*******************************************************************|************************************************************//
void DiagnosticsFile::tag_(const std::string &name,
                    const uint &column,
                    const std::string &statistic,
                    const std::string &system,
                    const uint &components)
{
  
  file_ << "<field column=\"" << column
        << "\" name=\"" << name
        << "\" statistic=\"" << statistic << "\"";
  if(!system.empty())                                                // is this part of a system?
  {
    file_ << " system=\"" << system << "\"";
  }
  if (components > 0)                                                // does it have subcomponents (i.e. is it rank>0)? 
  {
    file_ << " components=\"" << components << "\"";
  }
  file_ << " />" << std::endl << std::flush;

  
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to timestepping
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_timestep_()
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  file_ << (*bucket_).timestep_count() << " ";  
  file_ << (*bucket_).current_time() << " ";
  file_ << (*bucket_).elapsed_walltime() << " ";
  file_ << (*bucket_).timestep() << " ";
  
  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// check if the file_ is open
//*******************************************************************|************************************************************//
const bool DiagnosticsFile::is_open() const
{
  return file_.is_open();
}

//*******************************************************************|************************************************************//
// close the file_ (if open)
//*******************************************************************|************************************************************//
void DiagnosticsFile::close()
{
  if (is_open())
  {
    file_.close();                                                   // close the file_ member
  }
}

