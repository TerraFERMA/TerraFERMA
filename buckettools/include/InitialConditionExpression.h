
#ifndef __INITIAL_CONDITION_EXPRESSION_H
#define __INITIAL_CONDITION_EXPRESSION_H

#include <dolfin.h>
#include "BoostTypes.h"
#include "Bucket.h"

namespace buckettools
{

  // This class provides a method of looping over field initial conditions (as individual expressions)
  // and setting a mixed function space to those initial conditions
  class InitialConditionExpression : public dolfin::Expression
  {
  // only accessible within this class
  private:
    
    // a map from field component number to expression
    std::map< uint, Expression_ptr > expressions_;
  
  // accessible to everyone
  public:

    // Specific constructor (scalar)
    InitialConditionExpression(std::map< uint, Expression_ptr > expressions);
    
    // Specific constructor (vector)
    InitialConditionExpression(uint dim, std::map< uint, Expression_ptr > expressions);
    
    // no copy constructor for now
    
    // Default destructor (also calls DOLFIN expression destructor)
    virtual ~InitialConditionExpression();
    
    // overloaded DOLFIN expression eval
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const;
    
  };

}
#endif
