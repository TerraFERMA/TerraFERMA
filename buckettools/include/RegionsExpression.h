
#ifndef __REGIONS_EXPRESSION_H
#define __REGIONS_EXPRESSION_H

#include <dolfin.h>
#include "BoostTypes.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // RegionsExpression class:
  //
  // This class provides a method of overloading a dolfin Expression eval function by looping over cell ids and returning the results
  // of individual expressions in each of those regions.
  //*****************************************************************|************************************************************//
  class RegionsExpression : public dolfin::Expression
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//
    
    RegionsExpression(                                               // specific constructor (scalar)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const uint &dim,                               // specific constructor (vector)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const uint &dim0, const uint &dim1,            // specific constructor (tensor)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    RegionsExpression(const std::vector<size_t> &value_shape,          // specific constructor (alternate tensor)
                      std::map< std::size_t, Expression_ptr > expressions,
                      MeshFunction_size_t_ptr cell_ids);
    
    virtual ~RegionsExpression();                           // default destructor
    
    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    void eval(dolfin::Array<double>& values,                         // evaluate this expression at a point in a given cell
              const dolfin::Array<double>& x, 
              const ufc::cell &cell) const;
    
    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::map< std::size_t, Expression_ptr> expressions() const        // return the map of expressions (const version)
    { return expressions_; }
    
    std::map< std::size_t, Expression_ptr> expressions()                    // return the map of expressions (non-const version)
    { return expressions_; }
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available to this class
    
    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::size_t, Expression_ptr > expressions_;                   // map from region id to expression for that region
  
    MeshFunction_size_t_ptr cell_ids_;                                 // a (boost shared) pointer to a (uint) mesh function holding the cell ids

  };

}
#endif
