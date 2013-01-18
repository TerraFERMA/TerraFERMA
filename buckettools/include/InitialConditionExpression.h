
#ifndef __INITIAL_CONDITION_EXPRESSION_H
#define __INITIAL_CONDITION_EXPRESSION_H

#include <dolfin.h>
#include "BoostTypes.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // InitalConditionExpression class:
  //
  // This class provides a method of looping over field initial conditions (as individual expressions)
  // and setting a mixed function space to those initial conditions.
  // It is a derived dolfin expression and overloads much functionality in that base class
  //*****************************************************************|************************************************************//
  class InitialConditionExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    InitialConditionExpression(                                      // specific constructor (scalar)
                      std::map< std::size_t, Expression_ptr > expressions);
    
    InitialConditionExpression(const uint &dim,                      // specific constructor (vector)
                      std::map< std::size_t, Expression_ptr > expressions);
    
    
    virtual ~InitialConditionExpression();                           // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate this expression at a point in a given cell
              const dolfin::Array<double>& x, 
              const ufc::cell &cell) const;
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::size_t, Expression_ptr > expressions_;                   // map from component to initial condition expression for a function
  
  };

}
#endif
