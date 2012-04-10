
#include "StatisticsFile.h"
#include "Bucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
StatisticsFile::StatisticsFile(const std::string &name) : DiagnosticsFile(name)
{
                                                                     // do nothing... all handled by DiagnosticsFile constructor
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
StatisticsFile::~StatisticsFile()
{
                                                                     // do nothing... all handled by DiagnosticsFile destructor
}

//*******************************************************************|************************************************************//
// write a header for the model described in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::write_header(const Bucket &bucket)
{
  bucket.copy_diagnostics(bucket_);

  uint column = 1;                                                   // keep count of how many columns there are
  
  file_ << "<header>" << std::endl;                                  // initialize header xml
  header_constants_();                                               // write constant tags
  header_timestep_(column);                                          // write tags for the timesteps
  header_bucket_(column);                                            // write tags for the actual bucket variables - fields etc.
  file_ << "</header>" << std::endl << std::flush;                   // finalize header xml
}

//*******************************************************************|************************************************************//
// write data for the model described in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::write_data()
{
  
  data_timestep_();                                                 // write the timestepping information
  data_bucket_();                                                   // write the bucket data
  
  file_ << std::endl << std::flush;                                 // flush the buffer
  
}

//*******************************************************************|************************************************************//
// write a header for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::header_bucket_(uint &column)
{

  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin();        // loop over the systems
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
//    header_system_((*sys_it).second, column);                        // write the header for the system itself

    header_func_((*(*sys_it).second).fields_begin(),                // write the header for the fields in the system
                          (*(*sys_it).second).fields_end(), column);

    header_func_((*(*sys_it).second).coeffs_begin(),                // write the header for the coefficients in the system
                          (*(*sys_it).second).coeffs_end(), column);
  }

}

//*******************************************************************|************************************************************//
// write a header for the model systems in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::header_system_(const SystemBucket_ptr sys_ptr, 
                                     uint &column)
{
  tag_((*sys_ptr).name(), column++, "max");
  tag_((*sys_ptr).name(), column++, "min");
}

//*******************************************************************|************************************************************//
// write a header for a set of model fields
//*******************************************************************|************************************************************//
void StatisticsFile::header_func_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end, 
                                  uint &column)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given functions
                                                      f_it++)
  {
    if ((*(*f_it).second).include_in_statistics())                   // check they should be included
    {                                                                // yes, then populate header with default stats (min and max)
      if ((*(*(*f_it).second).function()).value_rank()==0)           // scalar (no components)
      {
        tag_((*(*f_it).second).name(), column++, "max", 
                              (*(*(*f_it).second).system()).name());
        tag_((*(*f_it).second).name(), column++, "min", 
                              (*(*(*f_it).second).system()).name());
        if ((*(*f_it).second).residualfunction())                    // all fields should get in here
        {
          tag_((*(*f_it).second).name(), column++, "res_max", 
                                (*(*(*f_it).second).system()).name());
          tag_((*(*f_it).second).name(), column++, "res_min", 
                                (*(*(*f_it).second).system()).name());
        }
      }
      else if ((*(*(*f_it).second).function()).value_rank()==1)      // vector (value_size components)
      {
        int components = (*(*(*f_it).second).function()).value_size();
        tag_((*(*f_it).second).name(), column, "max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        if ((*(*f_it).second).residualfunction())                    // all fields should get in here
        {
          tag_((*(*f_it).second).name(), column, "res_max", 
                    (*(*(*f_it).second).system()).name(), components);
          column+=components;
          tag_((*(*f_it).second).name(), column, "res_min", 
                    (*(*(*f_it).second).system()).name(), components);
          column+=components;
        }
      }
      else if ((*(*(*f_it).second).function()).value_rank()==2)      // tensor (value_dimension product components)
      {
        int components = 
          (*(*(*f_it).second).function()).value_dimension(0)*(*(*(*f_it).second).function()).value_dimension(1);
        tag_((*(*f_it).second).name(), column, "max", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        tag_((*(*f_it).second).name(), column, "min", 
                  (*(*(*f_it).second).system()).name(), components);
        column+=components;
        if ((*(*f_it).second).residualfunction())                    // all fields should get in here
        {
          tag_((*(*f_it).second).name(), column, "res_max", 
                    (*(*(*f_it).second).system()).name(), components);
          column+=components;
          tag_((*(*f_it).second).name(), column, "res_min", 
                    (*(*(*f_it).second).system()).name(), components);
          column+=components;
        }
      }
      else                                                           // unknown rank
      {
        dolfin::error("In StatisticsFile::header_bucket_, unknown function rank.");
      }

      header_functional_((*f_it).second, 
                        (*(*f_it).second).functionals_begin(),       // write header for any functionals associated with this field
                        (*(*f_it).second).functionals_end(), column);
    }
  }
}

//*******************************************************************|************************************************************//
// write a header for a set of model functionals
//*******************************************************************|************************************************************//
void StatisticsFile::header_functional_(const FunctionBucket_ptr f_ptr, 
                                        Form_const_it f_begin, 
                                        Form_const_it f_end, 
                                        uint &column)
{
  for (Form_const_it f_it = f_begin; f_it != f_end; f_it++)          // loop over the functional forms associated with the given
  {                                                                  // function bucket
    tag_((*f_ptr).name(), column++, (*f_it).first,                   // write tags for each functional
                                  (*(*f_ptr).system()).name());
  }
}

//*******************************************************************|************************************************************//
// write data for the model systems, fields and coefficients in the given bucket
//*******************************************************************|************************************************************//
void StatisticsFile::data_bucket_()
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  for (SystemBucket_const_it sys_it = (*bucket_).systems_begin(); 
                          sys_it != (*bucket_).systems_end(); sys_it++)
  {
//    data_system_((*sys_it).second);

    data_field_((*(*sys_it).second).fields_begin(), 
                                (*(*sys_it).second).fields_end());

    data_coeff_((*(*sys_it).second).coeffs_begin(), 
                                (*(*sys_it).second).coeffs_end());
  }

  file_.unsetf(std::ios::scientific);
  
}

//*******************************************************************|************************************************************//
// write data for a system
//*******************************************************************|************************************************************//
void StatisticsFile::data_system_(const SystemBucket_ptr sys_ptr)
{
  file_ << (*(*(*sys_ptr).function()).vector()).max() << " ";
  file_ << (*(*(*sys_ptr).function()).vector()).min() << " ";
}

//*******************************************************************|************************************************************//
// write data for a set of fields
//*******************************************************************|************************************************************//
void StatisticsFile::data_field_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given fields
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_statistics())                   // check if they should be included in the diagnostics
    {                                                                // yes, start with the default stats... min and max
      dolfin::Function func =                                        // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).function());
      dolfin::Function resfunc =                                     // take a deep copy of the subfunction so the vector is accessible
        *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).residualfunction());
      if (func.value_rank()==0)                                      // scalars (no components)
      {
        file_ << (*func.vector()).max() << " ";
        file_ << (*func.vector()).min() << " ";

        file_ << (*resfunc.vector()).max() << " ";
        file_ << (*resfunc.vector()).min() << " ";
      }
      else if (func.value_rank()==1)                                 // vectors (multiple components)
      {
        int components = func.value_size();
        for (uint i = 0; i < components; i++)
        {
          dolfin::Function funccomp = func[i];                                        // take a deep copy of the component of the subfunction
          file_ << (*funccomp.vector()).max() << " ";                // maximum for all components
        }

        for (uint i = 0; i < components; i++)
        {
          dolfin::Function funccomp = func[i];                       // take a deep copy of the component of the subfunction
          file_ << (*funccomp.vector()).min() << " ";                // minimum for all components
        }

        for (uint i = 0; i < components; i++)
        {
          dolfin::Function resfunccomp = resfunc[i];
          file_ << (*resfunccomp.vector()).max() << " ";             // maximum for all components
        }

        for (uint i = 0; i < components; i++)
        {
          dolfin::Function resfunccomp = resfunc[i];                 // take a deep copy of the component of the subfunction
          file_ << (*resfunccomp.vector()).min() << " ";             // minimum for all components
        }
      }
      else if (func.value_rank()==2)                                 // tensor (multiple components)
      {
        int dim0 = func.value_dimension(0);
        int dim1 = func.value_dimension(1);
        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function funccomp = func[i][j];                  // take a deep copy of the ijth component of the subfunction
            file_ << (*funccomp.vector()).max() << " ";              // maximum for all components
          }
        }

        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function funccomp = func[i][j];                  // take a deep copy of the ijth component of the subfunction
            file_ << (*funccomp.vector()).min() << " ";              // minimum for all components
          }
        }

        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function resfunccomp = resfunc[i][j];            // take a deep copy of the ijth component of the subfunction
            file_ << (*resfunccomp.vector()).max() << " ";           // maximum for all components
          }
        }

        for (uint i = 0; i < dim0; i++)
        {
          for (uint j = 0; j < dim1; j++)
          {
            dolfin::Function resfunccomp = resfunc[i][j];            // take a deep copy of the ijth component of the subfunction
            file_ << (*resfunccomp.vector()).min() << " ";           // minimum for all components
          }
        }
      }
      else                                                           // unknown rank
      {
        dolfin::error("In StatisticsFile::data_field_, unknown function rank.");
      }

      data_functional_((*(*f_it).second).functionals_begin(),        // wttie data for all functionals associated with this field
                            (*(*f_it).second).functionals_end());
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a set of coefficients
//*******************************************************************|************************************************************//
void StatisticsFile::data_coeff_(FunctionBucket_const_it f_begin, 
                                  FunctionBucket_const_it f_end)
{
  for (FunctionBucket_const_it f_it = f_begin; f_it != f_end;        // loop over the given coefficients
                                                            f_it++)
  {
    if ((*(*f_it).second).include_in_statistics())                  // check if this coefficient is to be included in diagnostics
    {
      if ((*(*f_it).second).type()=="Function")                      // this is a function coefficient so it has a vector member
      {
        dolfin::Function func =                                      // take a deep copy of the subfunction so the vector is accessible
          *boost::dynamic_pointer_cast< const dolfin::Function >((*(*f_it).second).function());
        if (func.value_rank()==0)                                    // scalars (no components)
        {
          file_ << (*func.vector()).max() << " ";
          file_ << (*func.vector()).min() << " ";
        }
        else if (func.value_rank()==1)                               // vectors (multiple components)
        {
          int components = func.value_size();
          for (uint i = 0; i < components; i++)
          {
            dolfin::Function funccomp = func[i];                     // take a deep copy of the component of the subfunction
            file_ << (*funccomp.vector()).max() << " ";              // maximum for all components
          }
          for (uint i = 0; i < components; i++)
          {
            dolfin::Function funccomp = func[i];                     // take a deep copy of the component of the subfunction
            file_ << (*funccomp.vector()).min() << " ";              // minimum for all components
          }
        }
        else if (func.value_rank()==2)                               // tensor (multiple components)
        {
          int dim0 = func.value_dimension(0);
          int dim1 = func.value_dimension(1);
          for (uint i = 0; i < dim0; i++)
          {
            for (uint j = 0; j < dim1; j++)
            {
              dolfin::Function funccomp = func[i][j];                // take a deep copy of the ijth component of the subfunction
              file_ << (*funccomp.vector()).max() << " ";            // maximum for all components
            }
          }
          for (uint i = 0; i < dim0; i++)
          {
            for (uint j = 0; j < dim1; j++)
            {
              dolfin::Function funccomp = func[i][j];                // take a deep copy of the ijth component of the subfunction
              file_ << (*funccomp.vector()).min() << " ";            // minimum for all components
            }
          }
        }
        else                                                         // unknown rank
        {
          dolfin::error("In StatisticsFile::data_coeff_, unknown function rank.");
        }
      }
      else                                                           // no vector available
      {
        Mesh_ptr mesh = (*(*(*f_it).second).system()).mesh();
        GenericFunction_ptr func = (*(*f_it).second).function();
        std::vector< double > values;
        (*func).compute_vertex_values(values, *mesh);
        if ((*func).value_rank()==0)                                 // scalars (no components)
        {
          file_ << *std::max_element(&values[0], &values[values.size()]) << " ";
          file_ << *std::min_element(&values[0], &values[values.size()]) << " ";
        }
        else if ((*func).value_rank()==1)                            // vectors (multiple components)
        {
          int components = (*func).value_size();
          for (uint i = 0; i < components; i++)
          {
            file_ << *std::max_element(&values[i*(*mesh).num_vertices()], 
                &values[(i+1)*(*mesh).num_vertices()]) << " ";// maximum for all components
          }
          for (uint i = 0; i < components; i++)
          {
            file_ << *std::min_element(&values[i*(*mesh).num_vertices()], 
                &values[(i+1)*(*mesh).num_vertices()]) << " ";// maximum for all components
          }
        }
        else if ((*func).value_rank()==2)                            // tensor (multiple components)
        {
          int dim0 = (*func).value_dimension(0);
          int dim1 = (*func).value_dimension(1);
          for (uint i = 0; i < dim0; i++)
          {
            for (uint j = 0; j < dim1; j++)
            {
              file_ << *std::max_element(&values[(2*j+i)*(*mesh).num_vertices()], 
                  &values[(2*j+i+1)*(*mesh).num_vertices()]) << " ";// maximum for all components
            }
          }
          for (uint i = 0; i < dim0; i++)
          {
            for (uint j = 0; j < dim1; j++)
            {
              file_ << *std::min_element(&values[(2*j+i)*(*mesh).num_vertices()], 
                  &values[(2*j+i+1)*(*mesh).num_vertices()]) << " ";// maximum for all components
            }
          }
        }
        else                                                         // unknown rank
        {
          dolfin::error("In StatisticsFile::data_coeff_, unknown function rank.");
        }
      }

      data_functional_((*(*f_it).second).functionals_begin(),        // write data for all functionals associated with this
                              (*(*f_it).second).functionals_end());  // coefficient
    }
  }
}

//*******************************************************************|************************************************************//
// write data for a set of functional forms
//*******************************************************************|************************************************************//
void StatisticsFile::data_functional_(Form_const_it s_begin, 
                                       Form_const_it s_end)
{
  for (Form_const_it s_it = s_begin; s_it != s_end; s_it++)          // loop over the given functionals
  {
    double statistic = dolfin::assemble((*(*s_it).second));          // assemble the functional
    file_ << statistic << " ";                                       // write to file
  }
}

