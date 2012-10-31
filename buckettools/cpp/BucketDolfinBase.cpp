
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

//*******************************************************************|************************************************************//
// add zeros to the diagonal on owned rows to ensure they remain in the sparsity pattern
//*******************************************************************|************************************************************//
void Assembler::add_zeros_diagonal(dolfin::GenericTensor& A)
{
  const uint rank = A.rank();
  assert(rank == 2);
  assert(A.size(0) == A.size(1));

  std::vector< const std::vector<uint>* > dofs(rank);
  std::vector< uint > global_row(1);
  for (uint i = 0; i < rank; ++i)
    dofs[i] = &global_row;
  const double zero = 0.0;

  const std::pair<uint, uint> row_range = A.local_range(0);
  const uint m = row_range.second - row_range.first;

  // Loop over rows
  for (uint row = 0; row < m; row++)
  {
    // Get global row number
    global_row[0] = row + row_range.first;
    A.add(&zero, dofs);

  }

}

bool abslessthan (const double &elem1, const double &elem2)
{
    return std::abs(elem1) < std::abs(elem2);
}


