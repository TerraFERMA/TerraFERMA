
#include "InitialConditionExpression.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

InitialConditionExpression::InitialConditionExpression(std::map< uint, GenericFunction_ptr > expressions) : dolfin::Expression(), expressions_(expressions)
{
  // Do nothing
}
InitialConditionExpression::InitialConditionExpression(uint dim, std::map< uint, GenericFunction_ptr > expressions) : dolfin::Expression(dim), expressions_(expressions)
{
  // Do nothing
}
InitialConditionExpression::~InitialConditionExpression()
{
  // should all be automatic from pyinst_ destructor
}

void InitialConditionExpression::eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
{
  double zero = 0.0;
  values=zero;
  for (std::map< uint, GenericFunction_ptr >::const_iterator expr = expressions_.begin();
                    expr != expressions_.end(); expr++)
  {
    dolfin::Array<double> tmp_values((*(*expr).second).value_size());
    (*(*expr).second).eval(tmp_values, x, cell);
    for (uint i = 0; i< tmp_values.size(); i++)
    {
      values[(*expr).first+i] = tmp_values[i];
    }
  }
  
}
