
#include "GenericDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
GenericDetectors::GenericDetectors() : number_detectors_(0), meshdim_(0), name_("uninitialized_string")
{
  // Do nothing
}

// Specific constructor
GenericDetectors::GenericDetectors(dolfin::uint number_detectors, dolfin::uint meshdim, std::string name) : number_detectors_(number_detectors), meshdim_(meshdim), name_(name)
{
  // Do nothing
}

// Default destructor
GenericDetectors::~GenericDetectors()
{
  clean_();
}
  
void GenericDetectors::clean_()
{
  while (!positions_.empty())
  {
    positions_.pop_back();
  }
}

void GenericDetectors::eval(std::vector< Array_double_ptr > &values,
                     dolfin::GenericFunction &function)
{
  Array_double_ptr value(new dolfin::Array<double>(function.value_size()));
  
  assert(values.empty());
  
  for (uint i = 0; i<positions_.size(); i++)
  {
    value.reset(new dolfin::Array<double>(function.value_size()));
    function.eval(*value, *positions_[i]);
    values.push_back(value);
  }
}

std::vector< Array_double_ptr >::iterator GenericDetectors::begin()
{
  return positions_.begin();
}

std::vector< Array_double_ptr >::const_iterator GenericDetectors::begin() const
{
  return positions_.begin();
}

std::vector< Array_double_ptr >::iterator GenericDetectors::end()
{
  return positions_.end();
}

std::vector< Array_double_ptr >::const_iterator GenericDetectors::end() const
{
  return positions_.end();
}

std::string GenericDetectors::name() const
{
  return name_;
}

uint GenericDetectors::dim() const
{
  return meshdim_;
}

uint GenericDetectors::size() const
{
  return number_detectors_;
}

std::string GenericDetectors::str() const
{
  std::stringstream s;
  
  for (dolfin::uint i = 0; i < positions_.size(); i++)
  {
    s << "detector " << i << std::endl;
    s << (*positions_[i]).str(true);
  }
  
  return s.str();
}
