
#include "Bucket.h"
#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

Bucket::Bucket()
{
  // Do nothing
}

Bucket::~Bucket()
{
  // Do nothing
}

void Bucket::register_detector(Detectors_ptr detector)
{
  detectors_.push_back(detector);
}

void Bucket::register_function(GenericFunction_ptr function)
{
  functions_.push_back(function);
}

std::vector< Detectors_ptr >::iterator Bucket::detectors_begin()
{
  return detectors_.begin();
}

std::vector< Detectors_ptr >::const_iterator Bucket::detectors_begin() const
{
  return detectors_.begin();
}

std::vector< Detectors_ptr >::iterator Bucket::detectors_end()
{
  return detectors_.end();
}

std::vector< Detectors_ptr >::const_iterator Bucket::detectors_end() const
{
  return detectors_.end();
}

std::vector< GenericFunction_ptr >::iterator Bucket::functions_begin()
{
  return functions_.begin();
}

std::vector< GenericFunction_ptr >::const_iterator Bucket::functions_begin() const
{
  return functions_.begin();
}

std::vector< GenericFunction_ptr >::iterator Bucket::functions_end()
{
  return functions_.end();
}

std::vector< GenericFunction_ptr >::const_iterator Bucket::functions_end() const
{
  return functions_.end();
}

void Bucket::clean_()
{
  while (!detectors_.empty())
  {
    detectors_.pop_back();
  }
  while (!functions_.empty())
  {
    functions_.pop_back();
  }
}
