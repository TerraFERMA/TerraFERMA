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


#include "DetectorsFile.h"
#include "Bucket.h"
#include "MPIBase.h"
#include "Logger.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
DetectorsFile::DetectorsFile(const std::string name, 
                             const MPI_Comm &comm, 
                             const Bucket *bucket) : DiagnosticsFile(name, comm, bucket)
{
  if (dolfin::MPI::size(mpicomm_)>1)
  {
#ifdef HAS_MPI                                                       // presumably true as size returned > 1
    int mpierr;
    mpierr = MPI_File_open(mpicomm_, (char*)(name_+".dat").c_str(),  // delete the previous dat file - FIXME: ugly hack
                           MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, 
                           MPI_INFO_NULL, &mpifile_);
    mpi_err(mpierr);
    mpierr = MPI_File_close(&mpifile_);
    mpi_err(mpierr);

    mpierr = MPI_File_open(mpicomm_, (char*)(name_+".dat").c_str(), 
                           MPI_MODE_CREATE + MPI_MODE_RDWR, 
                           MPI_INFO_NULL, &mpifile_);
    if (mpierr!=MPI_SUCCESS)
    {
      tf_err("MPI error opening MPI_File.", "MPI error: %d", mpierr);
    }
#endif
  }

  mpiwritecount_ = 0;                                                // incremented at every data dump
#ifdef HAS_MPI
  mpiwritelocation_ = 0;
#endif

}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
DetectorsFile::~DetectorsFile()
{
  if (dolfin::MPI::size(mpicomm_)>1)
  {
#ifdef HAS_MPI                                                       // presumably true as size return > 1
    int mpierr = MPI_File_close(&mpifile_);
    if (mpierr!=MPI_SUCCESS)
    {
      tf_err("MPI error closing MPI_File.", "MPI error: %d", mpierr);
    }
#endif
  }
}

//*******************************************************************|************************************************************//
// write header for the model described in the given bucket
//*******************************************************************|************************************************************//
void DetectorsFile::write_header()
{
  header_open_();
  header_constants_(dolfin::MPI::size(mpicomm_)>1);
  header_timestep_();
  header_bucket_();
  header_close_();

}

//*******************************************************************|************************************************************//
// write data for the model described in the attached bucket
//*******************************************************************|************************************************************//
void DetectorsFile::write_data()
{
  
  data_timestep_();
  data_bucket_();
  
  if (dolfin::MPI::size(mpicomm_)>1)
  {
    mpiwritecount_++;
    // quick sanity check...
#ifdef HAS_MPI
    int mpierr;
    MPI_Aint lb, doublesize;
    mpierr = MPI_Type_get_extent(MPI_DOUBLE_PRECISION, &lb, &doublesize);
    mpi_err(mpierr);
    assert(mpiwritelocation_==mpiwritecount_*ncolumns_*doublesize);
#endif
  }
  else
  {
    data_endlineflush_();
  }
  
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to timestepping
//*******************************************************************|************************************************************//
void DetectorsFile::data_timestep_()
{
  
  if (dolfin::MPI::size(mpicomm_)==1)
  {
    DiagnosticsFile::data_timestep_();
  }
  else
  {
#ifdef HAS_MPI
    int mpierr;
    double value;

    MPI_Aint lb, doublesize;
    mpierr = MPI_Type_get_extent(MPI_DOUBLE_PRECISION, &lb, &doublesize);
    mpi_err(mpierr);

    value = (double)(*bucket_).timestep_count();                     // we recast this to make it easier to count columns
    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      mpierr = MPI_File_write_at(mpifile_, mpiwritelocation_, 
                                 &value, 1, 
                                 MPI_DOUBLE_PRECISION, 
                                 MPI_STATUS_IGNORE);
      mpi_err(mpierr);
    }
    mpiwritelocation_ += doublesize;

    value = (*bucket_).current_time();
    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      mpierr = MPI_File_write_at(mpifile_, mpiwritelocation_, 
                                 &value, 1, 
                                 MPI_DOUBLE_PRECISION, 
                                 MPI_STATUS_IGNORE);
      mpi_err(mpierr);
    }
    mpiwritelocation_ += doublesize;

    value = (*bucket_).elapsed_walltime();
    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      mpierr = MPI_File_write_at(mpifile_, mpiwritelocation_, 
                                 &value, 1, 
                                 MPI_DOUBLE_PRECISION, 
                                 MPI_STATUS_IGNORE);
      mpi_err(mpierr);
    }
    mpiwritelocation_ += doublesize;

    value = (*bucket_).timestep();
    if (dolfin::MPI::rank(mpicomm_)==0)
    {
      mpierr = MPI_File_write_at(mpifile_, mpiwritelocation_, 
                                 &value, 1, 
                                 MPI_DOUBLE_PRECISION, 
                                 MPI_STATUS_IGNORE);
      mpi_err(mpierr);
    }
    mpiwritelocation_ += doublesize;
#endif
  }
  
}


//*******************************************************************|************************************************************//
// write header for the model described in the given bucket
//*******************************************************************|************************************************************//
void DetectorsFile::header_bucket_()
{
  for ( GenericDetectors_const_it d_it = (*bucket_).detectors_begin(); 
                                  d_it != (*bucket_).detectors_end(); 
                                  d_it++ )                           // loop over the detectors
  {
    header_detector_((*d_it).second);
  }

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();    // loop over the systems
                             sys_it != (*bucket_).systems_end(); 
                             sys_it++)
  {

    for ( FunctionBucket_const_it f_it = (*(*sys_it).second).fields_begin();// loop over the fields
                                  f_it != (*(*sys_it).second).fields_end(); 
                                  f_it++)
    {
      header_func_((*f_it).second);                                  // write the header for the fields in the system
    }


    for ( FunctionBucket_const_it f_it = (*(*sys_it).second).coeffs_begin();// loop over the fields
                                  f_it != (*(*sys_it).second).coeffs_end(); 
                                  f_it++)
    {
      header_func_((*f_it).second);                                  // write the header for the fields in the system
    }

  }
  
}

//*******************************************************************|************************************************************//
// write header for the detectors
//*******************************************************************|************************************************************//
void DetectorsFile::header_detector_(const GenericDetectors_ptr d_ptr)
{
  std::stringstream buffer;

  detectors_.push_back(d_ptr);

  for (uint dim = 0; dim<(*d_ptr).dim(); dim++)
  {
    buffer.str("");
    buffer << "position_" << dim;                                  // describe the detector positions
    tag_((*d_ptr).name(), buffer.str(), 
                                    "", (*d_ptr).size());
  }
  
}

//*******************************************************************|************************************************************//
// write header for the interaction between the detectors and the given functions
//*******************************************************************|************************************************************//
void DetectorsFile::header_func_(const FunctionBucket_ptr f_ptr)
{
  std::stringstream buffer;

  if ((*f_ptr).include_in_detectors())
  {

    functions_.push_back(f_ptr);

    for (std::vector<GenericDetectors_ptr>::const_iterator 
                                 d_it = detectors_.begin(); 
                                 d_it != detectors_.end(); 
                                 d_it++)
    {
      if ((*f_ptr).rank()==0)
      {
        tag_((*f_ptr).name(), (*(*d_it)).name(), 
              (*(*f_ptr).system()).name(), (*(*d_it)).size());
      }
      else if ((*f_ptr).rank()==1)
      {
        for (uint dim = 0; dim<(*f_ptr).size(); dim++)
        {
          buffer.str("");
          buffer << (*f_ptr).name() << "_" << dim;
          tag_(buffer.str(), (*(*d_it)).name(), 
                (*(*f_ptr).system()).name(), (*(*d_it)).size());
        }
      }
      else if ((*f_ptr).rank()==2)
      {
        for (uint dim0 = 0; dim0<(*f_ptr).dimension(0); dim0++)
        {
          for (uint dim1 = 0; dim1<(*f_ptr).dimension(1); dim1++)
          {
            buffer.str("");
            buffer << (*f_ptr).name() << "_" << dim0 << "_" << dim1;
            tag_(buffer.str(), (*(*d_it)).name(), 
                 (*(*f_ptr).system()).name(), (*(*d_it)).size());
          }
        }
      }
      else
      {
        tf_err("Unknown function rank.", "Rank: %d", (*f_ptr).rank());
      }
    }
  }
  
}

//*******************************************************************|************************************************************//
// write data for the attached bucket
//*******************************************************************|************************************************************//
void DetectorsFile::data_bucket_()
{
  
  for (std::vector<GenericDetectors_ptr>::const_iterator
                              d_it = detectors_.begin();
                              d_it != detectors_.end();
                              d_it++)
  {
    data_detector_(*d_it);
  }

  for (std::vector<FunctionBucket_ptr>::const_iterator
                              f_it = functions_.begin();
                              f_it != functions_.end();
                              f_it++)
  {
    data_func_(*f_it);
  }

}

//*******************************************************************|************************************************************//
// write data for the detectors
//*******************************************************************|************************************************************//
void DetectorsFile::data_detector_(const GenericDetectors_ptr d_ptr)
{

  const bool parallel = dolfin::MPI::size(mpicomm_)>1;
  const bool rank0 = dolfin::MPI::rank(mpicomm_)==0;                 // since all processors know the positions only rank 0 writes
                                                                     // perhaps all processes should write?
#ifdef HAS_MPI
  int mpierr;
  MPI_Aint lb, doublesize;
  mpierr = MPI_Type_get_extent(MPI_DOUBLE_PRECISION, &lb, &doublesize);
  mpi_err(mpierr);
#endif

  for (uint dim = 0; dim<(*d_ptr).dim(); dim++)
  {
    for (std::vector< Array_double_ptr >::const_iterator pos = 
                                    (*d_ptr).begin(); 
                          pos < (*d_ptr).end(); pos++)
    {   
      if (parallel)
      {
#ifdef HAS_MPI
        if(rank0)
        {
          mpierr = MPI_File_write_at(mpifile_, mpiwritelocation_, 
                                     &(**pos)[dim], 1, 
                                     MPI_DOUBLE_PRECISION, 
                                     MPI_STATUS_IGNORE);
          mpi_err(mpierr);
        }
        mpiwritelocation_ += doublesize;
#endif
      }
      else
      {
        data_((**pos)[dim]);
      }
    }
  }
  
}

//*******************************************************************|************************************************************//
// write data for the interaction between the detectors and the functions
//*******************************************************************|************************************************************//
void DetectorsFile::data_func_(const FunctionBucket_ptr f_ptr)
{
  
  const bool parallel = dolfin::MPI::size(mpicomm_)>1;
  
#ifdef HAS_MPI
  int mpierr;
  MPI_Aint lb, doublesize;
  mpierr = MPI_Type_get_extent(MPI_DOUBLE_PRECISION, &lb, &doublesize);
  mpi_err(mpierr);
#endif

  for ( std::vector<GenericDetectors_ptr>::const_iterator
                                         d_it = detectors_.begin();
                                         d_it != detectors_.end(); 
                                         d_it++)
  {
    std::vector< Array_double_ptr > values;
    
    GenericFunction_ptr func = (*f_ptr).iteratedfunction();

    (*(*d_it)).eval(values, *func, (*(*f_ptr).system()).mesh());
    std::vector< int > ids = (*(*d_it)).detector_ids((*(*f_ptr).system()).mesh());
    assert(values.size()==ids.size());
    
    for (uint dim = 0; dim < (*func).value_size(); dim++)
    {
      for(uint i=0; i < values.size(); i++)
      {
        if (parallel)
        {
#ifdef HAS_MPI
          MPI_Offset location = mpiwritelocation_ 
                              + (dim*((*(*d_it)).size()) 
                                 + ids[i])*doublesize;
          mpierr = MPI_File_write_at(mpifile_, location, 
                                     &(*values[i])[dim], 1, 
                                     MPI_DOUBLE_PRECISION, 
                                     MPI_STATUS_IGNORE);
          mpi_err(mpierr);
#endif
        }
        else
        {
          data_((*values[i])[dim]);
        }
      }
    }

    if (parallel)
    {
#ifdef HAS_MPI
      mpiwritelocation_ += (*func).value_size()*(*(*d_it)).size()*doublesize;
#endif
    }

  }
  
}

