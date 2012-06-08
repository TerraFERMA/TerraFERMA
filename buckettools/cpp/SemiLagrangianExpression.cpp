#include "SemiLagrangianExpression.h"
#include "BoostTypes.h"
#include "Bucket.h"
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
SemiLagrangianExpression::SemiLagrangianExpression(const Bucket *bucket, 
                                                   const double_ptr time,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &function,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &outside) :
                                                      dolfin::Expression(), 
                                                      bucket_(bucket), 
                                                      time_(time), 
                                                      funcname_(function),
                                                      velname_(velocity),
                                                      outname_(outside),
                                                      initialized_(false)
{
                                                                     // do nothing
}
    
//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
SemiLagrangianExpression::SemiLagrangianExpression(const uint &dim,
                                                   const Bucket *bucket, 
                                                   const double_ptr time,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &function,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &outside) :
                                                      dolfin::Expression(dim), 
                                                      bucket_(bucket), 
                                                      time_(time), 
                                                      funcname_(function),
                                                      velname_(velocity),
                                                      outname_(outside),
                                                      initialized_(false)
{
                                                                     // do nothing
}
    
//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
SemiLagrangianExpression::SemiLagrangianExpression(const uint &dim0, const uint &dim1,
                                                   const Bucket *bucket, 
                                                   const double_ptr time,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &function,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &outside) :
                                                      dolfin::Expression(dim0, dim1), 
                                                      bucket_(bucket), 
                                                      time_(time), 
                                                      funcname_(function),
                                                      velname_(velocity),
                                                      outname_(outside),
                                                      initialized_(false)
{
                                                                     // do nothing
}
    
//*******************************************************************|************************************************************//
// specific constructor (alternate tensor)
//*******************************************************************|************************************************************//
SemiLagrangianExpression::SemiLagrangianExpression(const std::vector<uint> &value_shape,
                                                   const Bucket *bucket, 
                                                   const double_ptr time,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &function,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                                                   const std::pair< std::string, std::pair< std::string, std::string > > &outside) :
                                                      dolfin::Expression(value_shape), 
                                                      bucket_(bucket), 
                                                      time_(time), 
                                                      funcname_(function),
                                                      velname_(velocity),
                                                      outname_(outside),
                                                      initialized_(false)
{
                                                                     // do nothing
}
    
//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SemiLagrangianExpression::~SemiLagrangianExpression()
{
  delete ufccellstar_;
  ufccellstar_ = NULL;
  delete[] xstar_;
  xstar_ = NULL;
  delete v_;
  v_ = NULL;
  delete oldv_;
  oldv_ = NULL;
  delete vstar_;
  vstar_ = NULL;
}

//*******************************************************************|************************************************************//
// initialize the expression
//*******************************************************************|************************************************************//
void SemiLagrangianExpression::init()
{
  if (!initialized_)
  {
    initialized_ = true;
    
    if(funcname_.second.first=="field")
    {
      func_ = (*(*(*bucket()).fetch_system(funcname_.first)).
               fetch_field(funcname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }
    else
    {
      func_ = (*(*(*bucket()).fetch_system(funcname_.first)).
               fetch_coeff(funcname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }

    if(velname_.second.first=="field")
    {
      vel_ = (*(*(*bucket()).fetch_system(velname_.first)).
               fetch_field(velname_.second.second)).
               genericfunction_ptr((*bucket()).current_time_ptr());
      oldvel_ = (*(*(*bucket()).fetch_system(velname_.first)).
               fetch_field(velname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }
    else
    {
      vel_ = (*(*(*bucket()).fetch_system(velname_.first)).
               fetch_coeff(velname_.second.second)).
               genericfunction_ptr((*bucket()).current_time_ptr());
      oldvel_ = (*(*(*bucket()).fetch_system(velname_.first)).
               fetch_coeff(velname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }

    if(outname_.second.first=="field")
    {
      out_ = (*(*(*bucket()).fetch_system(outname_.first)).
               fetch_field(outname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }
    else
    {
      out_ = (*(*(*bucket()).fetch_system(outname_.first)).
               fetch_coeff(outname_.second.second)).
               genericfunction_ptr((*bucket()).old_time_ptr());
    }

    assert(value_rank()==(*func_).value_rank());
    assert(value_rank()==(*out_).value_rank());
    for (uint i = 0; i < value_rank(); i++)
    {
      assert(value_dimension(i)==(*func_).value_dimension(i));
      assert(value_dimension(i)==(*out_).value_dimension(i));
    }

    mesh_ = (*(*bucket()).fetch_system(funcname_.first)).mesh();// use the mesh from the function system
    dim_ = (*mesh_).geometry().dim();

    assert((*vel_).value_rank()<2);
    if((*vel_).value_rank()==0)
    {
      assert(dim_==1);
    }
    else
    {
      assert((*vel_).value_dimension(0)==dim_);
    }

    ufccellstar_  = new dolfin::UFCCell(*mesh_);
    dolfincellit_ = new dolfin::CellIterator(*mesh_);
    xstar_        = new double[dim_];
    v_            = new dolfin::Array<double>(dim_);
    oldv_         = new dolfin::Array<double>(dim_);
    vstar_        = new dolfin::Array<double>(dim_);

  }
}
  
//*******************************************************************|************************************************************//
// overload dolfin expression eval
//*******************************************************************|************************************************************//
void SemiLagrangianExpression::eval(dolfin::Array<double>& values, 
                                    const dolfin::Array<double>& x, 
                                    const ufc::cell &cell) const
{
      
  bool outside = findpoint_(x, cell);
  
  dolfin::Array<double> xstar(dim_, xstar_);
  if (outside)
  {
    (*out_).eval(values, xstar, *ufccellstar_);
  }
  else
  {
    (*func_).eval(values, xstar, *ufccellstar_);
  }
      
}
    
//*******************************************************************|************************************************************//
// find the launch point
//*******************************************************************|************************************************************//
const bool SemiLagrangianExpression::findpoint_(const dolfin::Array<double>& x, 
                                                const ufc::cell &cell) const
{
  int cell_index;

  findvstar_(x, cell);

  for (uint k = 0; k < 2; k++)
  {
    for (uint i = 0; i < dim_; i++)
    {
      xstar_[i] = x[i] - ((*bucket()).timestep()/2.0)*(*vstar_)[i];
    }

    dolfin::Point p(dim_, xstar_);
    cell_index = (*mesh_).intersected_cell(p);
    if (cell_index<0)
    {
      return true;
    }
    else
    {
      (*ufccellstar_).update((*dolfincellit_)[cell_index]);
      dolfin::Array<double> xstar(dim_, xstar_);
      findvstar_(xstar, *ufccellstar_);
    }
  }

  for (uint i = 0; i < dim_; i++)
  {
    xstar_[i] = x[i] - ((*bucket()).timestep())*(*vstar_)[i];
  }

  dolfin::Point p(dim_, xstar_);
  cell_index = (*mesh_).intersected_cell(p);
  if (cell_index<0)
  {
    return true;
  }
  else
  {
    (*ufccellstar_).update((*dolfincellit_)[cell_index]);
  }
  return false;
}

void SemiLagrangianExpression::findvstar_(const dolfin::Array<double>& x, 
                                                const ufc::cell &cell) const
{
  (*vel_).eval(*v_, x, cell);
  (*oldvel_).eval(*oldv_, x, cell);
  for (uint i = 0; i < dim_; i++)
  {
    (*vstar_)[i] = 0.5*( (*v_)[i] + (*oldv_)[i] );
  }
}
