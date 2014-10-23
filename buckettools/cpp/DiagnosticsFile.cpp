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
#include "MPIBase.h"
#include "Usage.h"
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
DiagnosticsFile::DiagnosticsFile(const std::string &name, 
                                 const MPI_Comm &comm,
                                 const Bucket *bucket) : 
                                 name_(name), mpicomm_(comm), bucket_(bucket), ncolumns_(0)
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.open((char*)name.c_str());                                 // open the file_ member
  }
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
DiagnosticsFile::~DiagnosticsFile()
{
  close();                                                           // close the file_ member
}

//*******************************************************************|************************************************************//
// write opening lines of the xml header
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_open_()
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "<header>" << std::endl;                                // initialize header xml
  }
}

//*******************************************************************|************************************************************//
// write closing lines of the xml header
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_close_()
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "</header>" << std::endl << std::flush;                 // finalize header xml
  }
}

//*******************************************************************|************************************************************//
// write lines of the xml header for constants that do not vary throughout a simulation
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_constants_(const bool &binary)
{

  std::stringstream buffer;
  const int maxlencbuffer=1024;
  char cbuffer[maxlencbuffer];
  int cerr;

  constant_tag_("GitHash", "string", githash());                 // the git sha
  
  buffer.str("");
  buffer << __DATE__;
  buffer << " ";
  buffer << __TIME__;
  constant_tag_("CompileTime", "string", buffer.str());            // the compilation time
  
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

  if (binary)
  {
    constant_tag_("format", "string", "binary");

    int doublesize, intsize;
#ifdef HAS_MPI
    // currently binary output is only used for parallel detectors so the mpi definition should be used
    int mpierr;
    MPI_Aint mpidoublesize, mpiintsize;
    mpierr = MPI_Type_extent(MPI_DOUBLE_PRECISION, &mpidoublesize);
    mpi_err(mpierr);
    doublesize = (int)mpidoublesize;
    mpierr = MPI_Type_extent(MPI_INTEGER, &mpiintsize);
    mpi_err(mpierr);
    intsize = (int)mpiintsize;
#else
    doublesize = sizeof(double);
    intsize = sizeof(int);
#endif

    buffer.str(""); buffer << doublesize;
    constant_tag_("real_size", "integer", buffer.str());

    buffer.str(""); buffer << intsize;
    constant_tag_("integer_size", "integer", buffer.str());
  }
  else
  {
    constant_tag_("format", "string", "plain_text");
  }

}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to timestepping
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_timestep_()
{
  
  tag_("timestep", "value");                                         // the timestep count (i.e. number of timesteps taken)
  tag_("ElapsedTime", "value");                                      // the current time
  tag_("ElapsedWallTime", "value");                                  // the elapsed wall time
  tag_("dt", "value");                                               // the actual timestep
  
}

//*******************************************************************|************************************************************//
// write an xml tag for a constant that does not vary throughout a simulation
//*******************************************************************|************************************************************//
void DiagnosticsFile::constant_tag_(const std::string &name, 
                             const std::string &type, 
                             const std::string &value)
{
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "<constant name=\"" << name
          << "\" type=\"" << type
          << "\" value=\"" << value << "\" />" 
          << std::endl << std::flush;
  }

}

//*******************************************************************|************************************************************//
// write an xml tag for a variable in a simulation 
//*******************************************************************|************************************************************//
void DiagnosticsFile::tag_(const std::string &name,
                    const std::string &statistic,
                    const std::string &system,
                    const uint &components)
{
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << "<field column=\"" << ncolumns_+1
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

  if (components > 0)
  {
    ncolumns_+=components;
  }
  else
  {
    ncolumns_++;
  }
  
}

//*******************************************************************|************************************************************//
// end the line and flush
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_endlineflush_()
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << std::endl << std::flush;                               // flush the buffer
  }
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to timestepping
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_timestep_()
{
  
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.setf(std::ios::scientific);
    file_.precision(10);
    
    file_ << (*bucket_).timestep_count() << " ";  
    file_ << (*bucket_).current_time() << " ";
    file_ << (*bucket_).elapsed_walltime() << " ";
    file_ << (*bucket_).timestep() << " ";
    
    file_.unsetf(std::ios::scientific);
  }
  
}

//*******************************************************************|************************************************************//
// write generic data to file
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_(const int &value)
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_ << value << " ";
  }
}

//*******************************************************************|************************************************************//
// write generic data to file
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_(const double &value)
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.setf(std::ios::scientific);
    file_.precision(10);
    file_ << value << " ";
    file_.unsetf(std::ios::scientific);
  }
}

//*******************************************************************|************************************************************//
// write generic data to file
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_(const std::vector<double> &values)
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    file_.setf(std::ios::scientific);
    file_.precision(10);
    for (uint i = 0; i < values.size(); i++)
    {
      file_ << values[i] << " ";
    }
    file_.unsetf(std::ios::scientific);
  }
}

//*******************************************************************|************************************************************//
// close the file_ (if open)
//*******************************************************************|************************************************************//
void DiagnosticsFile::close()
{
  if (dolfin::MPI::rank(mpicomm_)==0)
  {
    if (file_.is_open())
    {
      file_.close();                                                 // close the file_ member
    }
  }
}

