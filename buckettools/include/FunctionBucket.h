
#ifndef __FUNCTIONBUCKET_H
#define __FUNCTIONBUCKET_H

#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  class SystemBucket;                                                // predeclare
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;        // so we can predeclare a pointer to it
  class FunctionBucket;                                              // predeclare the class itself
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // so we can predeclare a pointer to it
  
  enum solve_location { FUNCTIONBUCKET_FIELD, FUNCTIONBUCKET_COEFF };

  //*****************************************************************|************************************************************//
  // FunctionBucket class:
  //
  // The FunctionBucket class describes system functions and coefficients and provides data types
  // to the underlying functionals.
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

    const GenericFunction_ptr genericfunction_ptr(                   // return a constant (boost shared) pointer to the 
                                       const double_ptr time) const; // old or iterated function depending on the time pointer provided
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const double max(const double_ptr time, const int* index0, 
                                            const int* index1) const;// max of the function at the given time

    const double max(const double_ptr time, const int* index0) const
    { return max(time, index0, NULL); }

    const double max(const double_ptr time) const
    { return max(time, NULL, NULL); }

    const double min(const double_ptr time, const int* index0, 
                                            const int* index1) const;// min of the function at the given time

    const double min(const double_ptr time, const int* index0) const
    { return min(time, index0, NULL); }

    const double min(const double_ptr time) const
    { return min(time, NULL, NULL); }

    const double functionmax() const;                                // max of the function

    const double functionmin() const;                                // min of the function

    const std::string name() const                                   // return a constant string giving the function name
    { return name_; }

    const std::string type() const                                   // return a constant string giving the function type
    { return type_; }

    const std::string uflsymbol() const                              // return a constant string giving the ufl symbol 
    { return uflsymbol_; }                                           // for this function

    const uint index() const                                         // return a constant unsigned integer to the index of this
    { return index_; }                                               // function in the parent system

    SystemBucket* system()                                           // return a pointer to the parent system
    { return system_; }

    const SystemBucket* system() const                               // return a constant pointer to the parent system
    { return system_; }

    const FunctionSpace_ptr functionspace() const                    // return a constant (boost shared) pointer to the
    { return functionspace_; }                                       // functionspace
                                                                     // NOTE: if this is a field of a mixed system functionspace,
                                                                     // this will return a subspace


    const GenericFunction_ptr function() const                       // return a constant (boost shared) pointer to the 
    { return function_; }                                            // function 
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr oldfunction() const                    // return a constant (boost shared) pointer to the old
    { return oldfunction_; }                                         // function (previous timestep's values)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr iteratedfunction() const               // return a constant (boost shared) pointer to the iterated
    { return iteratedfunction_; }                                    // function (most up to date values within an iteration)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr changefunction() const                 // return a constant (boost shared) pointer to the change
    { return changefunction_; }                                      // function (change between timesteps)
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr residualfunction() const               // return a constant (boost shared) pointer to the residual
    { return residualfunction_; }                                    // function
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const GenericFunction_ptr snesupdatefunction() const             // return a constant (boost shared) pointer to the snes update
    { return snesupdatefunction_; }                                  // function
                                                                     // NOTE: if this is a field of a mixed
                                                                     // system functionspace, this will return a subfunction
                                                                     // so it will be necessary to make a deep copy to access
                                                                     // the vector

    const Expression_ptr icexpression() const                        // return a constant (boost shared) pointer to the initial
    { return icexpression_; }                                        // condition expression for this function

    const std::string change_normtype() const                        // return the change norm type
    { return change_normtype_; }

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    const double change();                                           // return the change (in the selected norm) in a field over a timestep

    const double functionalchange(Form_const_it f_it);               // return the change in the value of the given functional

    const double oldfunctionalvalue(Form_const_it f_it);             // return the old value of the given functional

    const double functionalvalue(Form_const_it f_it);                // return the value of the given functional

    void resetcalculated();                                          // reset the calculated booleans

    void refresh(const bool force=false);                            // refresh this function bucket - this may call solvers so 
                                                                     // its not recommened to call loosely

    void update();                                                   // update the timelevels of the function if it is potentially time dependent

    void update_nonlinear();                                         // update this function if it is potentially nonlinear

    void update_timedependent();                                     // update the function if it is potentially time dependent

    void cap_values();                                               // cap the values in the system vector associated with a field

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_functional_coeffs();                                 // attach the coefficients to the functionals of this function

    virtual void copy_diagnostics(FunctionBucket_ptr &function, 
                                  SystemBucket_ptr &system) const;   // copy the data necessary for the diagnostics data file(s)

    //***************************************************************|***********************************************************//
    // Functional data access
    //***************************************************************|***********************************************************//

    void register_functional(Form_ptr functional, 
                                           const std::string &name); // register a functional with the given name in this function

    Form_ptr fetch_functional(const std::string &name);              // return a (boost shared) pointer to a functional with the given name

    Form_it functionals_begin();                                     // return an iterator to the beginning of the functionals of this function

    Form_const_it functionals_begin() const;                         // return a constant iterator to the beginning of the functionals of this function

    Form_it functionals_end();                                       // return an iterator to the end of the functionals of this function

    Form_const_it functionals_end() const;                           // return a constant iterator to the end of the functionals of this function

    //***************************************************************|***********************************************************//
    // BC data access
    //***************************************************************|***********************************************************//

    void register_bcexpression(Expression_ptr bcexpression,          // register an expression for a bc in this function
                                          const std::string &name);

    Expression_ptr fetch_bcexpression(const std::string &name);      // return a (boost shared) pointer to a bc expression with the given name

    void register_bc(BoundaryCondition_ptr bc, 
                                        const std::string &name);    // register a bc in this function

    BoundaryCondition_it bcs_begin();                                // return an iterator to the beginning of the bcs of this function

    BoundaryCondition_const_it bcs_begin() const;                    // return a constant iterator to the beginning of the bcs of this function

    BoundaryCondition_it bcs_end();                                  // return an iterator to the end of the bcs of this function

    BoundaryCondition_const_it bcs_end() const;                      // return a constant iterator to the end of the bcs of this function

    int_BoundaryCondition_it orderedbcs_begin();                     // return an iterator to the beginning of the ordered bcs of this function

    int_BoundaryCondition_const_it orderedbcs_begin() const;         // return a constant iterator to the beginning of the ordered bcs of this function

    int_BoundaryCondition_it orderedbcs_end();                       // return an iterator to the end of the ordered bcs of this function

    int_BoundaryCondition_const_it orderedbcs_end() const;           // return a constant iterator to the end of the ordered bcs of this function

    void register_dirichletbc(BoundaryCondition_ptr bc, 
                                        const std::string &name);    // register a Dirichlet bc in this function

    BoundaryCondition_it dirichletbcs_begin();                       // return an iterator to the beginning of the bcs of this function

    BoundaryCondition_const_it dirichletbcs_begin() const;           // return a constant iterator to the beginning of the bcs of this function

    BoundaryCondition_it dirichletbcs_end();                         // return an iterator to the end of the bcs of this function

    BoundaryCondition_const_it dirichletbcs_end() const;             // return a constant iterator to the end of the bcs of this function

    int_BoundaryCondition_it ordereddirichletbcs_begin();            // return an iterator to the beginning of the ordered bcs of this function

    int_BoundaryCondition_const_it ordereddirichletbcs_begin() const;// return a constant iterator to the beginning of the ordered bcs of this function

    int_BoundaryCondition_it ordereddirichletbcs_end();              // return an iterator to the end of the ordered bcs of this function

    int_BoundaryCondition_const_it ordereddirichletbcs_end() const;  // return a constant iterator to the end of the ordered bcs of this function

    //***************************************************************|***********************************************************//
    // Reference point data access
    //***************************************************************|***********************************************************//

    void register_point(ReferencePoints_ptr point, 
                                        const std::string &name);    // register a point in this function

    ReferencePoints_it points_begin();                                // return an iterator to the beginning of the points of this function

    ReferencePoints_const_it points_begin() const;                    // return a constant iterator to the beginning of the points of this function

    ReferencePoints_it points_end();                                  // return an iterator to the end of the points of this function

    ReferencePoints_const_it points_end() const;                      // return a constant iterator to the end of the points of this function

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    void output(const bool &write_vis);                              // output diagnostics about this function

    virtual const bool include_in_visualization() const;             // return a boolean indicating if this function is included in 
                                                                     // visualization output

    virtual const bool include_residual_in_visualization() const;    // return a boolean indicating if the residual of this function is included in 
                                                                     // visualization output

    virtual const bool include_in_statistics() const;                // return a boolean indicating if this function is included in 
                                                                     // diagnostic output

    virtual const bool include_in_steadystate() const;               // return a boolean indicating if this function is included in 
                                                                     // steady state output

    virtual const bool include_functional_in_steadystate(const std::string &name) const;// return a boolean indicating if the named functional is included in 
                                                                     // steady state output

    virtual const bool include_in_detectors() const;                 // return a boolean indicating if this function is included in 
                                                                     // detectors output

    virtual const std::string str() const                            // return a string describing the contents of this function
    { return str(0); }

    virtual const std::string str(int indent) const;                 // return an indented string describing the contents 
                                                                     // of this function

    virtual const std::string functionals_str() const                // return a string describing the functionals of this function
    { return functionals_str(0); }

    virtual const std::string functionals_str(int indent) const;     // return an indented string describing the functionals 
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

    GenericFunction_ptr function_, oldfunction_, iteratedfunction_;  // (boost shared) pointers to different timelevel values of this function

    GenericFunction_ptr changefunction_;                             // (boost shared) pointer to the change in the function over a timestep

    GenericFunction_ptr residualfunction_;                           // (boost shared) pointer to the residual in the function (only fields)

    GenericFunction_ptr snesupdatefunction_;                         // (boost shared) pointer to the snes update of the function (only fields that use snes monitors)

    Expression_ptr icexpression_;                                    // (boost shared) pointer to an expression describing the initial condition

    int size_;                                                       // size of the function (most relevant for rank 1, vectors)

    std::vector< int > shape_;                                       // shape of the function (most relevant for rank 2, tensors)

    std::string rank_;                                               // a *string* describing the rank of the function
    
    int functiontype_;                                               // a *integer* describing the type of function bucket (field or coefficient)

    std::string type_;                                               // a *string* describing the type of function (function, expression, constant)

    Expression_ptr coefficientfunction_;                             // an expression used to set the values of a coefficient functional

    Form_ptr constantfunctional_;                                    // a functional that can be used to set a constant function
    
    File_ptr pvdfile_, respvdfile_;                                  // (boost shared) pointer to a pvd file output for this function

    double_ptr change_;                                              // change in the function in a norm

    bool_ptr change_calculated_;                                     // indicate if the change has been recalculated recently

    std::string change_normtype_;                                    // norm type to evaluate the change in

    double_ptr lower_cap_, upper_cap_;                               // caps on a field

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, Form_ptr > functionals_;                  // map from functional names to form (boost shared) pointers

    std::map< std::string, double_ptr > functional_values_,          // map from functional names to functional values
                                        oldfunctional_values_;

    std::map< std::string, bool_ptr > functional_calculated_;        // map from functional names to a boolean indicating if the
                                                                     // functional's been evaluated this timestep

    std::map< std::string, Expression_ptr > bcexpressions_;          // map from bc names to bc expression (boost shared) pointers
    
    std::map< std::string, BoundaryCondition_ptr > bcs_;             // map from bc::id names to (boost shared) pointers to bcs
    
    std::map< int, BoundaryCondition_ptr > orderedbcs_;              // map from int to (boost shared) pointers to bcs

    std::map< std::string, BoundaryCondition_ptr > dirichletbcs_;    // map from bc::id names to (boost shared) pointers to dirichlet bcs
    
    std::map< int, BoundaryCondition_ptr > ordereddirichletbcs_;     // map from int to (boost shared) pointers to dirichlet bcs

    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    virtual void checkpoint_options_();                              // checkpoint the options system for the functionbucket

    std::map< std::string, ReferencePoints_ptr > points_;            // map from bc::id names to (boost shared) pointers to bcs
    
    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                                   // empty the data structures of this function bucket

  };

  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;    // define a (boost shared) pointer to the function bucket class type

}
#endif
