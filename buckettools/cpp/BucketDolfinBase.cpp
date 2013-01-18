
#include "BucketDolfinBase.h"
#include <dolfin.h>

using namespace buckettools;

Side::Side(const uint &component, const double &side) : component_(component), side_(side)
{
  // Do nothing
}

Side::~Side()
{
  // Do nothing
}

bool Side::inside(const dolfin::Array<double>& x, bool on_boundary) const
{
  return (std::fabs(x[component_] - side_) < DOLFIN_EPS && on_boundary);
}

bool abslessthan (const double &elem1, const double &elem2)
{
    return std::abs(elem1) < std::abs(elem2);
}


