
#include "Detectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

Detectors::Detectors() : number_detectors_(0), meshdim_(0), name_("uninitialized_string")
{
  // Do nothing
}

Detectors::Detectors(dolfin::uint number_detectors, dolfin::uint meshdim, std::string name) : number_detectors_(number_detectors), meshdim_(meshdim), name_(name)
{
  // Do nothing
}

Detectors::~Detectors()
{
  clean_();
}
  
void Detectors::clean_()
{
  while (!positions_.empty())
  {
    positions_.pop_back();
  }
}

void Detectors::eval(std::vector< Array_double_ptr > &values,
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

std::vector< Array_double_ptr >::iterator Detectors::begin()
{
  return positions_.begin();
}

std::vector< Array_double_ptr >::const_iterator Detectors::begin() const
{
  return positions_.begin();
}

std::vector< Array_double_ptr >::iterator Detectors::end()
{
  return positions_.end();
}

std::vector< Array_double_ptr >::const_iterator Detectors::end() const
{
  return positions_.end();
}

std::string Detectors::name() const
{
  return name_;
}

uint Detectors::dim() const
{
  return meshdim_;
}

uint Detectors::size() const
{
  return number_detectors_;
}

std::string Detectors::str() const
{
  std::stringstream s;
  
  for (dolfin::uint i = 0; i < positions_.size(); i++)
  {
    s << "detector " << i << std::endl;
    s << (*positions_[i]).str(true);
  }
  
  return s.str();
}
