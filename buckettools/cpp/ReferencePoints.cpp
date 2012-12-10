
#include "ReferencePoints.h"
#include <dolfin.h>
#include <string>
#include <limits>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
ReferencePoints::ReferencePoints() : name_("uninitialized_string")
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
ReferencePoints::ReferencePoints(const Array_double_ptr coord, const FunctionSpace_ptr functionspace, 
                                          const std::string &name) : name_(name), functionspace_(functionspace)
{
  init_(coord);                                                      // initialize the detector 
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
ReferencePoints::ReferencePoints(const std::vector<double> &coord, const FunctionSpace_ptr functionspace, 
                                          const std::string &name) : name_(name), functionspace_(functionspace)
{
  init_(coord);                                                      // initialize the detector
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
ReferencePoints::~ReferencePoints()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// apply to a matrix
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericMatrix& A) const
{
  apply(&A, 0, 0);
}

//*******************************************************************|************************************************************//
// apply to the rhs of a linear problem
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericVector& b) const
{
  apply(0, &b, 0);
}

//*******************************************************************|************************************************************//
// apply to the matrix and rhs of a linear problem
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericMatrix& A, dolfin::GenericVector& b) const
{
  apply(&A, &b, 0);
}

//*******************************************************************|************************************************************//
// apply to the rhs of a nonlinear problem
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericVector& b, const dolfin::GenericVector& x) const
{
  apply(0, &b, &x);
}

//*******************************************************************|************************************************************//
// apply to the matrix and rhs of a nonlinear problem
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericMatrix& A,
                            dolfin::GenericVector& b,
                            const dolfin::GenericVector& x) const
{
  apply(&A, &b, &x);
}

//*******************************************************************|************************************************************//
// apply to the matrix and rhs of a nonlinear problem (implementation)
//*******************************************************************|************************************************************//
void ReferencePoints::apply(dolfin::GenericMatrix* A,
                            dolfin::GenericVector* b,
                            const dolfin::GenericVector* x) const
{
  if (dof_.size()==0)
  {
    return;                                                          // nothing to do
  }

  check_arguments_(A, b, x);                                         // check arguments

  const double spring = std::numeric_limits<double>::max()*std::numeric_limits<double>::epsilon();

  const uint size = dof_.size();
  std::vector<double> values(size, 0.0);

  if (x)                                                             // deal with nonlinear problems
  {
    (*x).get_local(&values[0], dof_.size(), &dof_[0]);
  }

  dolfin::log(dolfin::PROGRESS, "Applying reference points to linear system.");

  if (b)
  {
    if (A)
    {
      for (uint i = 0; i < size; i++)
      {
        values[i] = values[i]*spring;
      }
    }

    (*b).set(&values[0], size, &dof_[0]);
    (*b).apply("insert");
  }

  if (A)
  {
    std::vector< const std::vector<dolfin::DolfinIndex>* > block_dofs(2);
    std::vector< dolfin::DolfinIndex > row(1);
    for (uint i = 0; i < 2; ++i )
    {
      block_dofs[i] = &row;
    }
    
    for (uint i = 0; i < size; i++)
    {
      row[0] = dof_[i];
      (*A).add(&spring, block_dofs);
    }
    
    (*A).apply("add");
  }
}

//*******************************************************************|************************************************************//
// check the validity of the arguments to apply
//*******************************************************************|************************************************************//
void ReferencePoints::check_arguments_(dolfin::GenericMatrix* A, dolfin::GenericVector* b,
                                        const dolfin::GenericVector* x) const
{
  assert(functionspace_);

                                                                     // Check matrix and vector dimensions
  if (A && x && (*A).size(0) != (*x).size())
  {
    dolfin::error("Matrix dimension (%d rows) does not match vector dimension (%d) for application of reference points",
                 (*A).size(0), (*x).size());
  }

  if (A && b && (*A).size(0) != (*b).size())
  {
    dolfin::error("Matrix dimension (%d rows) does not match vector dimension (%d) for application of reference points",
                 (*A).size(0), (*b).size());
  }

  if (x && b && (*x).size() != (*b).size())
  {
    dolfin::error("Vector dimension (%d rows) does not match vector dimension (%d) for application of reference points",
                 (*x).size(), (*b).size());
  }

                                                                     // Check dimension of function space
  if (A && (*A).size(0) < (*functionspace_).dim())
  {
    dolfin::error("Dimension of function space (%d) too large for application of reference points to linear system (%d rows)",
                 (*functionspace_).dim(), (*A).size(0));
  }

  if (x && (*x).size() < (*functionspace_).dim())
  {
    dolfin::error("Dimension of function space (%d) too large for application of reference points linear system (%d rows)",
                 (*functionspace_).dim(), (*x).size());
  }

  if (b && (*b).size() < (*functionspace_).dim())
  {
    dolfin::error("Dimension of function space (%d) too large for application of reference points linear system (%d rows)",
                 (*functionspace_).dim(), (*b).size());
  }

}

//*******************************************************************|************************************************************//
// intialize a reference point from a (boost shared) pointer to a dolfin array
//*******************************************************************|************************************************************//
void ReferencePoints::init_(const Array_double_ptr coord)
{
  if(position_)
  {
    dolfin::error("In ReferencePoints::init_ intializing already initialized reference point.");
  }
  
  position_ = coord;

  const dolfin::Mesh& mesh = *(*functionspace_).mesh();
  const dolfin::GenericDofMap& dofmap = *(*functionspace_).dofmap();

  const std::size_t gdim = mesh.geometry().dim();

  double* pos = (*position_).data();
  uint dim = (*position_).size();
  assert(dim==gdim);
  const dolfin::Point point(dim, pos);
  int cellid = mesh.intersected_cell(point);

  if (cellid==-1)
  {
    if (dolfin::MPI::num_processes() == 1)
    {
      dolfin::log(dolfin::WARNING, "Failed to find cell at requested reference point coordinates.");
    }
  }
  else
  {
    const dolfin::Cell cell(mesh, cellid);

    boost::multi_array<double, 2> coordinates(boost::extents[dofmap.cell_dimension(cellid)][gdim]);
    dofmap.tabulate_coordinates(coordinates, cell);

    const std::vector<dolfin::DolfinIndex>& cell_dofs = dofmap.cell_dofs(cellid);

    std::vector<double> dist(dofmap.cell_dimension(cellid), 0.0);

    for (uint i = 0; i < dofmap.cell_dimension(cellid); ++i)
    {
      for (uint j = 0; j < gdim; ++j)
      {
        dist[i] += std::pow((coordinates[i][j] - (*position_)[j]), 2);
      }
    }

    double mindist = *std::min_element(&dist[0], &dist[dist.size()]);
    for (uint i = 0; i < dofmap.cell_dimension(cellid); ++i)
    {
      if (std::fabs(dist[i]-mindist) < DOLFIN_EPS)
      {
        dof_.push_back(cell_dofs[i]);
      }
    }
    
  }

}

//*******************************************************************|************************************************************//
// intialize a reference point from a std vector
//*******************************************************************|************************************************************//
void ReferencePoints::init_(const std::vector<double> &coord)
{
  Array_double_ptr arraypoint(new dolfin::Array<double>(coord.size()));
  for (uint i = 0; i<coord.size(); i++)
  {
    (*arraypoint)[i] = coord[i];
  }
  
  init_(arraypoint);
}

