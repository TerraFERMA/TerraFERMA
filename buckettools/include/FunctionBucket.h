// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __FUNCTIONBUCKET_H
#define __FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include "PointDetectors.h"
#include "ReferencePoint.h"
#include <dolfin.h>

namespace buckettools
{

  class SystemBucket;                                                // predeclare
  typedef std::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  class FunctionBucket;                                              // predeclare the class itself
  typedef std::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // so we can predeclare a pointer to it
  
  enum function_type { FUNCTIONBUCKET_FIELD, FUNCTIONBUCKET_COEFF };

  //*****************************************************************|************************************************************//
  // FunctionBucket class:
  //
  // The FunctionBucket class describes system functions and coefficients.
  //*****************************************************************|************************************************************//
  class FunctionBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    FunctionBucket();                                                // default constructor

    FunctionBucket(SystemBucket* system);                            // specific constructor
    
    virtual ~FunctionBucket();                                       // default destructor

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const GenericFunction_ptr genericfunction_ptr(                   // return a constant (std shared) pointer to the 
                                       const double_ptr time) const; // old or iterated function depending on the time pointer provided
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr genericfunction_ptr(                   // similar to the other genericfunction_ptr subroutine 
                                 const std::string &function_type) const; // except here the function_type is more general and described
                                                                     // by a string

    void cachevector(const std::string &function_type);

    void clearcachedvector();

    const_PETScVector_ptr basevector(const std::string &function_type) const; // return a pointer to the base (potentially system) vector

    dolfin::PETScVector vector(const std::string &function_type,     // return a vector of values for this function bucket and
                                        const int &component) const; // the specified function_type

    dolfin::PETScVector vector(const std::string &function_type, 
                            const std::vector<int>* components=NULL) 
                                                              const;
    IS component_is(const int &component) const;

    IS components_is(const std::vector<int>* components=NULL) const;

    double max(const std::string &function_type, 
                     const uint component) const;

    double max(const std::string &function_type, 
                     const std::vector<int>* components=NULL) const;

    double min(const std::string &function_type, 
                     const uint component) const;

    double min(const std::string &function_type, 
                     const std::vector<int>* components=NULL) const;

    double norm(const std::string &function_type, 
                      const std::string &norm_type, 
                      const uint component) const;

    double norm(const std::string &function_type, 
                      const std::string &norm_type, 
                      const std::vector<int>* components=NULL) const;

    double change(const uint component);

    double change(const std::vector<int>* components=NULL);          // return the change (in the selected norm) in a field over a timestep

    const std::string name() const                                   // return a constant string giving the function name
    { return name_; }

    const std::string type() const                                   // return a constant string giving the function type
    { return type_; }

    const std::string uflsymbol() const                              // return a constant string giving the ufl symbol 
    { return uflsymbol_; }                                           // for this function

    const uint index() const                                         // return a constant unsigned integer to the index of this
    { return index_; }                                               // function in the parent system

    const std::size_t rank() const                                   // return the rank of the function
    { return shape_.size(); }

    const std::size_t size() const;                                  // return the size of the function function

    const std::vector< std::size_t > shape() const                   // return the shape of the function
    { return shape_; }

    const std::size_t dimension(const std::size_t &i) const;         // return the size of a given dimension

    const bool symmetric() const;                                    // return if this is a symmetric tensor field (false otherwise)

    SystemBucket* system()                                           // return a pointer to the parent system
    { return system_; }

    const SystemBucket* system() const                               // return a constant pointer to the parent system
    { return system_; }

    const FunctionSpace_ptr functionspace() const                    // return a constant (std shared) pointer to the
    { return functionspace_; }                                       // functionspace
                                                                     // NOTE: if this is a field of a mixed system functionspace,
                                                                     // this will return a subspace


    const GenericFunction_ptr function() const                       // return a constant (std shared) pointer to the 
    { return function_; }                                            // function 
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr oldfunction() const                    // return a constant (std shared) pointer to the old
    { return oldfunction_; }                                         // function (previous timestep's values)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr iteratedfunction() const               // return a constant (std shared) pointer to the iterated
    { return iteratedfunction_; }                                    // function (most up to date values within an iteration)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr changefunction() const                 // return a constant (std shared) pointer to the change
    { return changefunction_; }                                      // function (change between timesteps)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr residualfunction() const               // return a constant (std shared) pointer to the residual
    { return residualfunction_; }                                    // function
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr snesupdatefunction() const             // return a constant (std shared) pointer to the snes update
    { return snesupdatefunction_; }                                  // function
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const Expression_ptr icexpression() const                        // return a constant (std shared) pointer to the initial
    { return icexpression_; }                                        // condition expression for this function

    const std::string change_normtype() const                        // return the change norm type
    { return change_normtype_; }

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    void refresh(const bool force=false);                            // refresh this function bucket - this may call solvers so 
                                                                     // its not recommened to call loosely

    void update();                                                   // update the timelevels of the function if it is potentially time dependent

    void update_nonlinear();                                         // update this function if it is potentially nonlinear

    void update_timedependent();                                     // update the function if it is potentially time dependent

    void resetcalculated();                            

    void postprocess_values();                                       // cap the values in the system vector associated with a field

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_form_coeffs();                                       // attach the coefficients to the functionals of this function

    //***************************************************************|***********************************************************//
    // BC data access
    //***************************************************************|***********************************************************//

    void register_bcfunction(GenericFunction_ptr bcfunction,         // register a function for a bc in this function
                                          const std::string &name);

    GenericFunction_ptr fetch_bcfunction(const std::string &name);   // return a (std shared) pointer to a bc function with the given name

    void register_dirichletbc(DirichletBC_ptr bc, 
                                        const std::string &name);    // register a Dirichlet bc in this function

    DirichletBC_it dirichletbcs_begin();                             // return an iterator to the beginning of the bcs of this function

    DirichletBC_const_it dirichletbcs_begin() const;                 // return a constant iterator to the beginning of the bcs of this function

    DirichletBC_it dirichletbcs_end();                               // return an iterator to the end of the bcs of this function

    DirichletBC_const_it dirichletbcs_end() const;                   // return a constant iterator to the end of the bcs of this function

    //***************************************************************|***********************************************************//
    // Reference point data access
    //***************************************************************|***********************************************************//

    void register_referencepoint(ReferencePoint_ptr point, 
                                        const std::string &name);    // register a point in this function

    ReferencePoint_it referencepoints_begin();                       // return an iterator to the beginning of the points of this function

    ReferencePoint_const_it referencepoints_begin() const;           // return a constant iterator to the beginning of the points of this function

    ReferencePoint_it referencepoints_end();                         // return an iterator to the end of the points of this function

    ReferencePoint_const_it referencepoints_end() const;             // return a constant iterator to the end of the points of this function

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    virtual const bool include_in_visualization() const;             // return a boolean indicating if this function is included in 
                                                                     // visualization output

    virtual const bool include_residual_in_visualization() const;    // return a boolean indicating if the residual of this function is included in 
                                                                     // visualization output

    virtual const bool include_in_statistics() const;                // return a boolean indicating if this function is included in 
                                                                     // diagnostic output

    virtual const bool include_in_steadystate() const;               // return a boolean indicating if this function is included in 
                                                                     // steady state output

    virtual const bool include_in_detectors() const;                 // return a boolean indicating if this function is included in 
                                                                     // detectors output

    virtual const std::string str(int indent=0) const;               // return an indented string describing the contents 
                                                                     // of this function

    void checkpoint();                                               // checkpoint the functionbucket

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the function name

    std::string uflsymbol_;                                          // the function ufl symbol

    uint index_;                                                     // the index of this field or coefficient in the system
                                                                     // (most relevant for field which are subfunctions of the system function)

    SystemBucket* system_;                                           // a pointer to the parent system

    FunctionSpace_ptr functionspace_;                                // the functionspace (may be a subspace) of this function (if it is
                                                                     // a field or a coefficient function)

    GenericFunction_ptr function_, oldfunction_, iteratedfunction_;  // (std shared) pointers to different timelevel values of this function

    GenericFunction_ptr changefunction_;                             // (std shared) pointer to the change in the function over a timestep

    GenericFunction_ptr residualfunction_;                           // (std shared) pointer to the residual in the function (only fields)

    GenericFunction_ptr snesupdatefunction_;                         // (std shared) pointer to the snes update of the function (only fields that use snes monitors)

    Expression_ptr icexpression_;                                    // (std shared) pointer to an expression describing the initial condition

    std::vector< std::size_t > shape_;                               // shape of the function

    int functiontype_;                                               // a *integer* describing the type of function bucket (field or coefficient)

    std::string type_;                                               // a *string* describing the type of function (function, expression, constant)

    Expression_ptr coefficientfunction_;                             // an expression used to set the values of a coefficient function

    Form_ptr constantfunctional_;                                    // a functional that can be used to set a constant function
    
    double_ptr change_;                                              // change in the function in a norm

    bool_ptr change_calculated_;                                     // indicate if the change has been recalculated recently

    std::string change_normtype_;                                    // norm type to evaluate the change in

    std::vector<double_ptr> lower_cap_, upper_cap_;                  // caps on a field

    std::vector<IS> component_is_;                                   // IS for components (into system vector!)

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, GenericFunction_ptr > bcfunctions_;       // map from bc names to bc function (std shared) pointers
    
    ordered_map<const std::string, DirichletBC_ptr> dirichletbcs_;   // map from bc names to (std shared) pointers to dirichlet bcs

    ordered_map<const std::string, ReferencePoint_ptr> referencepoints_;   // map from bc::id names to (std shared) pointers to bcs

    std::vector< PointDetectors_ptr > zeropoints_;                   // a list of zero points

    const_PETScVector_ptr cachedvector_;                             // cache the values of the vector temporarily

    std::string cachedvectortype_;                                   // the cached vector type (if it exists)
    
    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill_is_();                                                 // fill the index sets for this function's components
 
    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    virtual void checkpoint_options_();                              // checkpoint the options system for the functionbucket

  };

  typedef std::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // define a (std shared) pointer to the function bucket class type

}
#endif
