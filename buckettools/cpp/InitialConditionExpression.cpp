
#include "InitialConditionExpression.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor (scalar)
//*******************************************************************|************************************************************//
InitialConditionExpression::InitialConditionExpression(std::map< uint, 
                                     Expression_ptr > expressions) : 
                                            dolfin::Expression(), 
                                            expressions_(expressions)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor (vector)
//*******************************************************************|************************************************************//
InitialConditionExpression::InitialConditionExpression(uint dim, std::map< uint, Expression_ptr > expressions) : dolfin::Expression(dim), expressions_(expressions)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
InitialConditionExpression::~InitialConditionExpression()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// overload dolfin expression eval to loop over a map (from component) of expressions and evaluate them all to a single expression 
// value array
//*******************************************************************|************************************************************//
void InitialConditionExpression::eval(dolfin::Array<double>& values, 
                                      const dolfin::Array<double>& x, 
                                      const ufc::cell &cell) const
{
  double zero = 0.0;                                                 // zero the full value array
  values=zero;
  for (std::map< uint, Expression_ptr >::const_iterator expr =       // loop over the expressions in the map
                    expressions_.begin();
                    expr != expressions_.end(); expr++)
  {
    dolfin::Array<double> tmp_values((*(*expr).second).value_size());// set up a temporary value array
    (*(*expr).second).eval(tmp_values, x, cell);                     // evaluate the component expression
    for (uint i = 0; i< tmp_values.size(); i++)                      // loop over the temporary value array
    {
      values[(*expr).first+i] = tmp_values[i];                       // insert component value into the full value array
    }
  }
  
}
