
#ifndef __INITIAL_CONDITION_EXPRESSION_H
#define __INITIAL_CONDITION_EXPRESSION_H

#include <dolfin.h>
#include "Bucket.h"

namespace buckettools
{

  class InitialConditionExpression : public dolfin::Expression
  {
  private:
    
    std::map< uint, GenericFunction_ptr > expressions_;
    
  public:
    InitialConditionExpression(std::map< uint, GenericFunction_ptr > expressions);
    
    InitialConditionExpression(uint dim, std::map< uint, GenericFunction_ptr > expressions);
    
    // no copy constructor for now
    
    virtual ~InitialConditionExpression();
    
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const;
    
  };

}
#endif