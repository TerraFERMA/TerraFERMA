
#include "GenericDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
GenericDetectors::GenericDetectors() : number_detectors_(0), meshdim_(0), name_("uninitialized_string")
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
GenericDetectors::GenericDetectors(const uint &number_detectors, 
                                   const uint &meshdim, 
                                   const std::string &name) : 
                                   number_detectors_(number_detectors), 
                                   meshdim_(meshdim), name_(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
GenericDetectors::~GenericDetectors()
{
  clean_();                                                          // empty the data structures
}
  
//*******************************************************************|************************************************************//
// base eval implementation - evaluates a function at the detector positions and returns a vector of arrays of values
//*******************************************************************|************************************************************//
void GenericDetectors::eval(std::vector< Array_double_ptr > &values,
                            const dolfin::GenericFunction &function,
                            Mesh_ptr mesh)
{
  Array_double_ptr 
            value(new dolfin::Array<double>(function.value_size())); // set up a pointer to be reset later to contain the
  
  assert(values.empty());                                            // check the values are empty

  double* pos;
  uint dim;
  int id;
  
  std::map< Mesh_ptr, std::vector< int > >::const_iterator c_it = cell_ids_.find(mesh);
  bool have_cells = ( c_it != cell_ids_.end() );
  std::vector< int > cellids;
  if (have_cells)
  {
    cellids = (*c_it).second;
  }

  for (uint i = 0; i<positions_.size(); i++)                         // loop over the detector positions
  {
    value.reset(new dolfin::Array<double>(function.value_size()));   // reset the value

    if( have_cells )
    {
      id = cellids[i];
    }
    else
    {
      pos = (*positions_[i]).data().get(); 
      dim = (*positions_[i]).size();
      const dolfin::Point point(dim, pos);
      id = (*mesh).intersected_cell(point);
      cellids.push_back(id);
    }

    if(id==-1)
    {
      dolfin::error("Unable to find cell in mesh.");
    }
    else
    {
      const dolfin::Cell cell(*mesh, id);
      const dolfin::UFCCell ufc_cell(cell);
      function.eval(*value, *positions_[i], ufc_cell);                 // use the dolfin eval to evaluate the function
      values.push_back(value);                                         // record the value
    }
  }

  if (!have_cells)
  {
    cell_ids_[mesh] = cellids;
  }

}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the positions_ vector
//*******************************************************************|************************************************************//
std::vector< Array_double_ptr >::iterator GenericDetectors::begin()
{
  return positions_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the positions_ vector
//*******************************************************************|************************************************************//
std::vector< Array_double_ptr >::const_iterator GenericDetectors::begin() const
{
  return positions_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the positions_ vector
//*******************************************************************|************************************************************//
std::vector< Array_double_ptr >::iterator GenericDetectors::end()
{
  return positions_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the positions_ vector
//*******************************************************************|************************************************************//
std::vector< Array_double_ptr >::const_iterator GenericDetectors::end() const
{
  return positions_.end();
}

//*******************************************************************|************************************************************//
// return a string describing the positions of the detectors
//*******************************************************************|************************************************************//
const std::string GenericDetectors::str() const
{
  std::stringstream s;
  
  for (uint i = 0; i < positions_.size(); i++)                       // loop over the positions
  {
    s << "detector " << i << std::endl;
    s << (*positions_[i]).str(true);                                 // use the dolfin array str output
  }
  
  return s.str();
}

//*******************************************************************|************************************************************//
// empty the data structures in the genericdetectors class
//*******************************************************************|************************************************************//
void GenericDetectors::clean_()
{
  while (!positions_.empty())                                        // empty the positions vector
  {
    positions_.pop_back();
  }
}

