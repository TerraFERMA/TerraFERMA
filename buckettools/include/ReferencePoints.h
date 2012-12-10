
#ifndef __REFERENCEPOINTS_H
#define __REFERENCEPOINTS_H

#include <dolfin.h>
#include "BoostTypes.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // ReferencePoints class:
  //
  // ReferencePoints implements a method for setting reference points on functions.
  //*****************************************************************|************************************************************//
  class ReferencePoints
  {
  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    ReferencePoints();                                                // default constructor
    
    ReferencePoints(const Array_double_ptr coord, const FunctionSpace_ptr functionspace,
                                    const std::string &name);        // specific constructor (dolfin array)
    
    ReferencePoints(const std::vector<double> &coord, const FunctionSpace_ptr functionspace,
                                        const std::string &name);    // specific constructor (std vector)
    
    ~ReferencePoints();                                              // default destructor
    
    //***************************************************************|***********************************************************//
    // Application
    //***************************************************************|***********************************************************//

    void apply(dolfin::GenericMatrix& A) const;                      // apply to a matrix

    void apply(dolfin::GenericVector& b) const;                      // apply to the rhs of a linear problem

    void apply(dolfin::GenericMatrix& A, 
                                     dolfin::GenericVector& b) const;// apply to the matrix and rhs of a linear problem

    void apply(dolfin::GenericVector& b, 
                               const dolfin::GenericVector& x) const;// apply to the rhs of a nonlinear problem

    void apply(dolfin::GenericMatrix& A,
               dolfin::GenericVector& b,
               const dolfin::GenericVector& x) const;                // apply to the matrix and rhs of a nonlinear problem

    void apply(dolfin::GenericMatrix* A,
               dolfin::GenericVector* b,
               const dolfin::GenericVector* x) const;                // apply to the matrix and rhs of a nonlinear problem (implementation)

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string      name_;                                          // reference point name (not really necessary)
    Array_double_ptr position_;                                      // array giving location of reference point
    FunctionSpace_ptr functionspace_;                                // functionspace to which this reference point is applied
    std::vector<std::size_t> dof_;                                  // list of the degrees of freedom 

    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_(const Array_double_ptr coord);                        // initialize the point detector at a point defined by an array
    
    void init_(const std::vector<double> &coord);                    // initialize the point detector at a point defined by a vector
    
    //***************************************************************|***********************************************************//
    // Application
    //***************************************************************|***********************************************************//

    void check_arguments_(dolfin::GenericMatrix* A, 
                          dolfin::GenericVector* b,
                          const dolfin::GenericVector* x) const;     // check the validity of the arguments to apply

    
  };

  typedef boost::shared_ptr< ReferencePoints > ReferencePoints_ptr;    // define a (boost shared) pointer for this class type

}

#endif
