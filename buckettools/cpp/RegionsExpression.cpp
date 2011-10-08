
#include "RegionsExpression.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(std::map< uint, 
                                     Expression_ptr > expressions,
                                     MeshFunction_uint_ptr cell_ids) : 
                                            dolfin::Expression(), 
                                            expressions_(expressions),
                                            cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const uint &dim, 
                                     std::map< uint, Expression_ptr > 
                                               expressions,
                                     MeshFunction_uint_ptr cell_ids) : 
                                     dolfin::Expression(dim), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (tensor)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const uint &dim0, const uint &dim1, 
                                     std::map< uint, Expression_ptr > 
                                               expressions,
                                     MeshFunction_uint_ptr cell_ids) : 
                                     dolfin::Expression(dim0, dim1), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (alternate tensor)
//*******************************************************************|************************************************************//
RegionsExpression::RegionsExpression(const std::vector<uint>
                                                        &value_shape, 
                                     std::map< uint, Expression_ptr > 
                                               expressions,
                                     MeshFunction_uint_ptr cell_ids) : 
                                     dolfin::Expression(value_shape), 
                                     expressions_(expressions),
                                     cell_ids_(cell_ids)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
RegionsExpression::~RegionsExpression()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// overload dolfin expression eval to loop over a map (from component) of expressions and evaluate them all to a single expression 
// value array
//*******************************************************************|************************************************************//
void RegionsExpression::eval(dolfin::Array<double>& values, 
                                      const dolfin::Array<double>& x, 
                                      const ufc::cell &cell) const
{
  uint id = (*cell_ids_)[cell.index];
  uint_Expression_const_it e_it = expressions_.find(id);
  if (e_it == expressions_.end())
  {
    dolfin::error("Unknown region id %d in RegionsExpression eval", id);
  }
  else
  {
    (*(*e_it).second).eval(values, x, cell);
  }
}

