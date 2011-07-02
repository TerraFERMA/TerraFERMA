
#include "PointDetectors.h"
#include "GenericDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

PointDetectors::PointDetectors() : GenericDetectors()
{
  // Do nothing
}

PointDetectors::PointDetectors(Array_double_ptr point, std::string name) : GenericDetectors(1, (*point).size(), name)
{
  init_(point);
}

PointDetectors::PointDetectors(std::vector<double> *point, 
                               std::string name) : GenericDetectors(1, (*point).size(), name)
{
  init_(point);
}

PointDetectors::~PointDetectors()
{
  // Do nothing - taken care of by base destructor
}

void PointDetectors::init_(Array_double_ptr point)
{
  if(!positions_.empty())
  {
    dolfin::error("In PointDetectors::init_ intializing already initialized detectors.");
  }
  
  positions_.push_back(point);
}

void PointDetectors::init_(std::vector<double> *point)
{
  if(!positions_.empty())
  {
    dolfin::error("In PointDetectors::init_ intializing already initialized detectors.");
  }
  
  Array_double_ptr arraypoint(new dolfin::Array<double>((*point).size()));
  for (uint i = 0; i<(*point).size(); i++)
  {
    (*arraypoint)[i] = (*point)[i];
  }
  
  positions_.push_back(arraypoint);
}
