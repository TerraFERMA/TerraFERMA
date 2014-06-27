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


#include "DiagnosticsFile.h"
#include "Bucket.h"
#include "MPIBase.h"
#include "builddefs.h"
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
                                 const MPI_Comm &comm) : 
                                 name_(name), mpicomm_(comm), ncolumns_(0)
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
// write lines of the xml header for constants that do not vary throughout a simulation
//*******************************************************************|************************************************************//
void DiagnosticsFile::header_constants_(const bool &binary)
{

  std::stringstream buffer;
  const int maxlencbuffer=1024;
  char cbuffer[maxlencbuffer];
  int cerr;

  constant_tag_("GitHash", "string", __GIT_SHA__);                 // the git sha
  
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
// write generic data for a function (field or coefficient function)
//*******************************************************************|************************************************************//
void DiagnosticsFile::data_function_(dolfin::Function &func, 
                                     const bool &include_norms)
{
  if (func.value_rank()==0)                                          // scalars (no components)
  {
    double max, min;
    double norml2, normlinf;

    max = (*func.vector()).max();
    min = (*func.vector()).min();
    if (include_norms)
    {
      norml2   = (*func.vector()).norm("l2");
      normlinf = (*func.vector()).norm("linf");
    }

    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      file_ << max << " ";
      file_ << min << " ";

      if (include_norms)
      {
        file_ << norml2 << " ";
        file_ << normlinf << " ";
      }
    }
  }
  else if (func.value_rank()==1)                                     // vectors (multiple components)
  {
    int components = func.value_size();

    std::vector<double> max(components, 0.0), min(components, 0.0);
    std::vector<double> norml2(components, 0.0), normlinf(components, 0.0);

    for (uint i = 0; i < components; i++)
    {
      dolfin::Function funccomp = func[i];                           // take a deep copy of the component of the subfunction
      max[i] = (*funccomp.vector()).max();
      min[i] = (*funccomp.vector()).min();
      if (include_norms)
      {
        norml2[i]   = (*funccomp.vector()).norm("l2");
        normlinf[i] = (*funccomp.vector()).norm("linf");
      }
    }

    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      for (uint i = 0; i < components; i++)
      {
        file_ << max[i] << " ";                                      // maximum for all components
      }
      for (uint i = 0; i < components; i++)
      {
        file_ << min[i] << " ";                                      // minimum for all components
      }
      if (include_norms)
      {
        for (uint i = 0; i < components; i++)
        {
          file_ << norml2[i] << " ";
        }
        for (uint i = 0; i < components; i++)
        {
          file_ << normlinf[i] << " ";
        }
      }
    }

  }
  else if (func.value_rank()==2)                                     // tensor (multiple components)
  {
    const bool symmetric = ((*(*func.function_space()).element()).num_sub_elements() != func.value_size());
    int dim0 = func.value_dimension(0);
    int dim1 = func.value_dimension(1);
    std::vector<double> max(dim0*dim1, 0.0), min(dim0*dim1, 0.0);
    std::vector<double> norml2(dim0*dim1, 0.0), normlinf(dim0*dim1, 0.0);
    for (uint i = 0; i < dim0; i++)
    {
      for (uint j = 0; j < dim1; j++)
      {
        std::size_t k; 
        if (symmetric) 
        {
          if (j >= i) 
          {
            k = i*dim1 + j - (i*(i+1))/2; 
          }
          else
          { 
            k = j*dim1 + i - (j*(j+1))/2; 
          }
        }
        else
        {
          k = i*dim1 + j;
        }
        dolfin::Function funccomp = func[k];                         // take a deep copy of the ijth component of the subfunction
        max[i*dim1 + j] = (*funccomp.vector()).max();
        min[i*dim1 + j] = (*funccomp.vector()).min();
        if (include_norms)
        {
          norml2[i*dim1 + j]   = (*funccomp.vector()).norm("l2");
          normlinf[i*dim1 + j] = (*funccomp.vector()).norm("linf");
        }
      }
    }
        
    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      for (uint i = 0; i < dim0*dim1; i++)
      {
        file_ << max[i] << " ";                                      // maximum for all components
      }
      for (uint i = 0; i < dim0*dim1; i++)
      {
        file_ << min[i] << " ";                                      // minimum for all components
      }
      if (include_norms)
      {
        for (uint i = 0; i < dim0*dim1; i++)
        {
          file_ << norml2[i] << " ";
        }
        for (uint i = 0; i < dim0*dim1; i++)
        {
          file_ << normlinf[i] << " ";
        }
      }
    }
  }
  else                                                               // unknown rank
  {
    dolfin::error("In StatisticsFile::data_function_, unknown function rank.");
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

