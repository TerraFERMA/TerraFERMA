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


#include "ConvergenceFile.h"
#include "Bucket.h"
#include "SystemBucket.h"
#include "SolverBucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
ConvergenceFile::ConvergenceFile(const std::string &name, 
                                 const std::string &systemname, 
                                 const std::string &solvername) : 
                                      DiagnosticsFile(name),
                                      systemname_(systemname),
                                      solvername_(solvername)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
ConvergenceFile::~ConvergenceFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  uint column = 1;                                                   // keep count of how many columns there are
  
  file_ << "<header>" << std::endl;                                  // initialize header xml
  header_constants_();                                               // write constant tags
  header_timestep_(column);                                          // write tags for the timesteps
  header_iteration_(column);                                         // write tags for the iterations
  header_bucket_(column);                                            // write tags for the actual bucket variables - fields etc.
  file_ << "</header>" << std::endl << std::flush;                   // finalize header xml
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::write_data()
{
  
  data_timestep_();                                                  // write the timestepping information
  data_iteration_();                                                 // write the iteration information
  data_bucket_();                                                    // write the bucket data
  
  file_ << std::endl << std::flush;                                  // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::header_iteration_(uint &column)
{
  
  tag_("NonlinearSystemsIteration", column++, "value");              // the nonlinear systems iteration
  tag_("NonlinearIteration", column++, "value");                     // the nonlinear solver iteration
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_bucket_(uint &column)
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);

  header_system_(sys_ptr, column);                                   // write the header for the system itself

  header_func_((*sys_ptr).fields_begin(),                            // write the header for the fields in the system
                        (*sys_ptr).fields_end(), column);

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::header_system_(const SystemBucket_ptr sys_ptr, 
                                     uint &column)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  tag_((*sys_ptr).name(), column++, "max");
  tag_((*sys_ptr).name(), column++, "min");
  tag_((*sys_ptr).name(), column++, "res_max");
  tag_((*sys_ptr).name(), column++, "res_min");
  tag_((*sys_ptr).name(), column++, "res_norm(l2)");
  tag_((*sys_ptr).name(), column++, "res_norm(inf)");
  if ((*sol_ptr).type()=="SNES")
  {
    tag_((*sys_ptr).name(), column++, "update_max");
    tag_((*sys_ptr).name(), column++, "update_min");
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void ConvergenceFile::header_func_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end, 
                                  uint &column)
{
  SolverBucket_ptr sol_ptr = (*(*(*f_begin).second).system()).fetch_solver(solvername_);

  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given functions
                                                      f_it++)
  {
    if ((*(*(*f_it).second).function()).value_rank()==0)             // scalar (no components)
    {
      tag_((*(*f_it).second).name(), column++, "max", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), column++, "min", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), column++, "res_max", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), column++, "res_min", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), column++, "res_norm(l2)", 
                            (*(*(*f_it).second).system()).name());
      tag_((*(*f_it).second).name(), column++, "res_norm(linf)", 
                            (*(*(*f_it).second).system()).name());
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), column++, "update_max", 
                              (*(*(*f_it).second).system()).name());
        tag_((*(*f_it).second).name(), column++, "update_min", 
                              (*(*(*f_it).second).system()).name());
      }
    }
    else if ((*(*(*f_it).second).function()).value_rank()==1)        // vector (value_size components)
    {
      int components = (*(*(*f_it).second).function()).value_size();
      tag_((*(*f_it).second).name(), column, "max", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "min", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_max", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_min", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_norm(l2)", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_norm(linf)", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), column, "update_max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "update_min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
    }
    else if ((*(*(*f_it).second).function()).value_rank()==2)        // tensor (value_dimension product components)
    {
      int components = 
        (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
      tag_((*(*f_it).second).name(), column, "max", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "min", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_max", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_min", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_norm(l2)", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      tag_((*(*f_it).second).name(), column, "res_norm(linf)", 
                (*(*(*f_it).second).system()).name(), components);
      column+=components;
      if ((*sol_ptr).type()=="SNES")
      {
        tag_((*(*f_it).second).name(), column, "update_max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "update_min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
      }
    }
    else                                                             // unknown rank
    {
      dolfin::error("In ConvergenceFile::header_bucket_, unknown function rank.");
    }

  }
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to iterations
//*******************************************************************|************************************************************//
void ConvergenceFile::data_iteration_()
{
  SolverBucket_ptr sol_ptr = (*(*bucket_).fetch_system(systemname_)).fetch_solver(solvername_);
  
  file_ << (*bucket_).iteration_count() << " ";  
  file_ << (*sol_ptr).iteration_count() << " ";
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void ConvergenceFile::data_bucket_()
{
  SystemBucket_ptr sys_ptr = (*bucket_).fetch_system(systemname_);
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  data_system_(sys_ptr);

  data_field_((*sys_ptr).fields_begin(), (*sys_ptr).fields_end());

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void ConvergenceFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  SolverBucket_ptr sol_ptr = (*sys_ptr).fetch_solver(solvername_);

  file_ << (*(*(*sys_ptr).iteratedfunction()).vector()).max() << " ";
  file_ << (*(*(*sys_ptr).iteratedfunction()).vector()).min() << " ";
  file_ << (*(*(*sys_ptr).residualfunction()).vector()).max() << " ";
  file_ << (*(*(*sys_ptr).residualfunction()).vector()).min() << " ";
  file_ << (*(*(*sys_ptr).residualfunction()).vector()).norm("l2") << " ";
  file_ << (*(*(*sys_ptr).residualfunction()).vector()).norm("linf") << " ";
  if ((*sol_ptr).type()=="SNES")
  {
    file_ << (*(*(*sys_ptr).snesupdatefunction()).vector()).max() << " ";
    file_ << (*(*(*sys_ptr).snesupdatefunction()).vector()).min() << " ";
  }
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void ConvergenceFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  SolverBucket_ptr sol_ptr = (*(*(*f_begin).second).system()).fetch_solver(solvername_);

  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    dolfin::Function func =                                          // take a deep copy of the subfunction so the vector is accessible
      *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).iteratedfunction());
    dolfin::Function resfunc =                                       // take a deep copy of the subfunction so the vector is accessible
      *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).residualfunction());
    if (func.value_rank()==0)                                        // scalars (no components)
    {
      file_ << (*func.vector()).max() << " ";
      file_ << (*func.vector()).min() << " ";

      file_ << (*resfunc.vector()).max() << " ";
      file_ << (*resfunc.vector()).min() << " ";
      file_ << (*resfunc.vector()).norm("l2") << " ";
      file_ << (*resfunc.vector()).norm("linf") << " ";

      if ((*sol_ptr).type()=="SNES")
      {
        dolfin::Function upfunc =                                   // take a deep copy of the subfunction so the vector is accessible
          *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).snesupdatefunction());

        file_ << (*upfunc.vector()).max() << " ";
        file_ << (*upfunc.vector()).min() << " ";
      }
    }
    else if (func.value_rank()==1)                                   // vectors (multiple components)
    {
      int components = func.value_size();
      for (uint i = 0; i < components; i++)
      {
        dolfin::Function funccomp = func[i];                         // take a deep copy of the component of the subfunction
        file_ << (*funccomp.vector()).max() << " ";                  // maximum for all components
      }

      for (uint i = 0; i < components; i++)
      {
        dolfin::Function funccomp = func[i];                         // take a deep copy of the component of the subfunction
        file_ << (*funccomp.vector()).min() << " ";                  // minimum for all components
      }

      for (uint i = 0; i < components; i++)
      {
        dolfin::Function resfunccomp = resfunc[i];
        file_ << (*resfunccomp.vector()).max() << " ";               // maximum for all components
      }

      for (uint i = 0; i < components; i++)
      {
        dolfin::Function resfunccomp = resfunc[i];                   // take a deep copy of the component of the subfunction
        file_ << (*resfunccomp.vector()).min() << " ";               // minimum for all components
      }

      for (uint i = 0; i < components; i++)
      {
        dolfin::Function resfunccomp = resfunc[i];
        file_ << (*resfunccomp.vector()).norm("l2") << " ";               // maximum for all components
      }

      for (uint i = 0; i < components; i++)
      {
        dolfin::Function resfunccomp = resfunc[i];                   // take a deep copy of the component of the subfunction
        file_ << (*resfunccomp.vector()).norm("linf") << " ";               // minimum for all components
      }

      if ((*sol_ptr).type()=="SNES")
      {
        dolfin::Function upfunc =                                   // take a deep copy of the subfunction so the vector is accessible
          *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).snesupdatefunction());

        for (uint i = 0; i < components; i++)
        {
          dolfin::Function upfunccomp = upfunc[i];
          file_ << (*upfunccomp.vector()).max() << " ";               // maximum for all components
        }

        for (uint i = 0; i < components; i++)
        {
          dolfin::Function upfunccomp = upfunc[i];                    // take a deep copy of the component of the subfunction
          file_ << (*upfunccomp.vector()).min() << " ";               // minimum for all components
        }
      }
    }
    else if (func.value_rank()==2)                                   // tensor (multiple components)
    {
      int dim0 = func.value_dimension(0);
      int dim1 = func.value_dimension(1);
      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function funccomp = func[i*dim1 + j];                    // take a deep copy of the ijth component of the subfunction
          file_ << (*funccomp.vector()).max() << " ";                // maximum for all components
        }
      }

      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function funccomp = func[i*dim1 + j];                    // take a deep copy of the ijth component of the subfunction
          file_ << (*funccomp.vector()).min() << " ";                // minimum for all components
        }
      }

      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function resfunccomp = resfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
          file_ << (*resfunccomp.vector()).max() << " ";             // maximum for all components
        }
      }

      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function resfunccomp = resfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
          file_ << (*resfunccomp.vector()).min() << " ";             // minimum for all components
        }
      }

      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function resfunccomp = resfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
          file_ << (*resfunccomp.vector()).norm("l2") << " ";             // maximum for all components
        }
      }

      for (uint i = 0; i < dim0; i++)
      {
        for (uint j = 0; j < dim1; j++)
        {
          dolfin::Function resfunccomp = resfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
          file_ << (*resfunccomp.vector()).norm("linf") << " ";             // minimum for all components
        }
      }

      if ((*sol_ptr).type()=="SNES")
      {
        dolfin::Function upfunc =                                    // take a deep copy of the subfunction so the vector is accessible
          *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).snesupdatefunction());

        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function upfunccomp = upfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
            file_ << (*upfunccomp.vector()).max() << " ";            // maximum for all components
          }
        }

        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function upfunccomp = upfunc[i*dim1 + j];              // take a deep copy of the ijth component of the subfunction
            file_ << (*upfunccomp.vector()).min() << " ";            // minimum for all components
          }
        }
      }
    }
    else                                                             // unknown rank
    {
      dolfin::error("In ConvergenceFile::data_field_, unknown function rank.");
    }
  }
}


