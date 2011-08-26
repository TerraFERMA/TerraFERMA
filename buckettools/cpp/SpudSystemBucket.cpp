
#include "BoostTypes.h"
#include "InitialConditionExpression.h"
#include "SpudSystemBucket.h"
#include "SystemsWrapper.h"
#include "SpudBase.h"
#include "SpudFunctionBucket.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudSystemBucket::SpudSystemBucket(std::string optionpath, 
                                            Bucket* bucket) : 
                                            optionpath_(optionpath), 
                                            SystemBucket(bucket)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudSystemBucket::~SpudSystemBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// fill the system bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSystemBucket::fill()
{
  base_fill_();                                                      // fill in the base data (could be called from the
                                                                     // constructor?)

  systemfunction_fill_();                                            // register the functionspace and system functions

  fields_fill_();                                                    // initialize the fields (subfunctions) of this system
 
  expcoeffs_fill_();                                                 // initialize the coefficient expressions (and constants)
                                                                     // (can't do coefficient functions now because it's unlikely we 
                                                                     // have all the coefficient functionspaces)

  solvers_fill_();                                                   // initialize the nonlinear solvers in this system

}

//*******************************************************************|************************************************************//
// loop over the coefficients and initialize any that are coefficient functions
//*******************************************************************|************************************************************//
void SpudSystemBucket::funccoeffs_fill()
{
  
  for (FunctionBucket_it f_it = coeffs_begin(); f_it != coeffs_end();// loop over all coefficients
                                                              f_it++) 
  {                                                                  // recast as a spud derived class and initialize
    (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_function_coeff();
  }                                                                  // (check that this is a coefficient function within this
                                                                     // function)

}

//*******************************************************************|************************************************************//
// return a string describing the contents of the spud system
//*******************************************************************|************************************************************//
const std::string SpudSystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << " (" 
                                << optionpath() << ")" << std::endl;
  indent++;
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// fill the system bucket base data 
//*******************************************************************|************************************************************//
void SpudSystemBucket::base_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // get the system name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/ufl_symbol";           // get the system ufl symbol
  serr = Spud::get_option(buffer.str(), uflsymbol_); 
  spud_err(buffer.str(), serr);

  std::string meshname;                                              // get the system mesh
  buffer.str(""); buffer << optionpath() << "/mesh/name";
  serr = Spud::get_option(buffer.str(), meshname); 
  spud_err(buffer.str(), serr);
  mesh_ = (*bucket_).fetch_mesh(meshname);                           // and extract it from the bucket

}

//*******************************************************************|************************************************************//
// fill the system functionspace and function data 
//*******************************************************************|************************************************************//
void SpudSystemBucket::systemfunction_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  functionspace_ = ufc_fetch_functionspace(name(), mesh());          // fetch the first functionspace we can grab from the ufc 
                                                                     // for this system

  function_.reset( new dolfin::Function(functionspace_) );           // declare the function on this functionspace
  buffer.str(""); buffer << name() << "::Function";
  (*function_).rename( buffer.str(), buffer.str() );

  oldfunction_.reset( new dolfin::Function(functionspace_) );        // declare the old function on this functionspace
  buffer.str(""); buffer << name() << "::OldFunction";
  (*oldfunction_).rename( buffer.str(), buffer.str() );

  iteratedfunction_.reset( new dolfin::Function(functionspace_) );   // declare the iterated function on this functionspace
  buffer.str(""); buffer << name() << "::IteratedFunction";
  (*iteratedfunction_).rename( buffer.str(), buffer.str() );

}

//*******************************************************************|************************************************************//
// fill in the data about each field (or subfunction) of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::fields_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer

                                                                     // prepare the system initial condition expression:
  uint component = 0;                                                // initialize a counter for the scalar components of this
                                                                     // system
  std::map< uint, Expression_ptr > icexpressions;                    // set up a map from component to initial condition expression

  buffer.str("");  buffer << optionpath() << "/field";               // find out how many fields we have
  int nfields = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfields; i++)                                 // loop over the fields in the options dictionary
  {
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]";

                                                                     // declare a new field function bucket assuming this system is
                                                                     // its parent
    SpudFunctionBucket_ptr field(new SpudFunctionBucket( buffer.str(), this ));
    (*field).field_fill(i);                                          // fill in this field (providing its index in the system)
    register_field(field, (*field).name());                          // register this field in the system bucket
                                  
                                                                     // insert the field's initial condition expression into a 
                                                                     // temporary system map:
    uint_Expression_it e_it = icexpressions.find(component);         // check if this component already exists
    if (e_it != icexpressions.end())
    {
      dolfin::error(                                                 // if it does, issue an error
      "IC Expression with component number %d already exists in icexpressions map.", 
                                                        component);
    }
    else
    {
      icexpressions[component] = (*field).icexpression();            // if it doesn't, insert it into the map
    }

    component += (*(*field).icexpression()).value_size();            // increment the component count by the size of this field
                                                                     // (i.e. no. of scalar components)
  }

  collect_bcs_();                                                    // collect all the bcs together for convenience later
  apply_ic_(component, icexpressions);                               // apply the initial condition to the system function
  apply_bc_();                                                       // apply the boundary conditions we just collected

}

//*******************************************************************|************************************************************//
// given a map from components to field initial condition expressions initialize the system with a combined initial condition
//*******************************************************************|************************************************************//
void SpudSystemBucket::apply_ic_(const uint &component, 
              const std::map< uint, Expression_ptr > &icexpressions)
{
  Expression_ptr ic;
  if (component==1)
  {
    ic.reset( new InitialConditionExpression(icexpressions) );       // the system function is scalar so set up a scalar ic expression
  }
  else
  {                                                                  // multiple components so set up a multi-component ic
                                                                     // expression
    ic.reset( new InitialConditionExpression(component, icexpressions));
  }
  (*oldfunction_).interpolate(*ic);                                  // interpolate the initial condition onto the old function
  (*iteratedfunction_).vector() = (*oldfunction_).vector();          // set the iterated function vector to the old function vector
  (*function_).vector() = (*oldfunction_).vector();                  // set the function vector to the old function vector
}

//*******************************************************************|************************************************************//
// apply the vector of system boundary conditions to the system function vectors to ensure consisten initial and boundary conditions
//*******************************************************************|************************************************************//
void SpudSystemBucket::apply_bc_()
{
  for (std::vector<BoundaryCondition_ptr>::const_iterator            // loop over the vector of bcs
                            bc = bcs_begin(); bc != bcs_end(); bc++)
  {
    (*(*bc)).apply((*oldfunction_).vector());
    (*(*bc)).apply((*iteratedfunction_).vector());
    (*(*bc)).apply((*function_).vector());
  }
}

//*******************************************************************|************************************************************//
// fill in the data about each coefficient expression (or constant) of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::expcoeffs_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  
  buffer.str("");  buffer << optionpath() << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());                    // find out how many coefficients we have
  for (uint i = 0; i < ncoeffs; i++)                                 // and loop over them
  {
    buffer.str(""); buffer << optionpath() << "/coefficient[" 
                                                        << i << "]";

                                                                     // initialize a new function bucket for this coefficient
                                                                     // (regardless of type!) assuming this system is its parent
    SpudFunctionBucket_ptr coeff( new SpudFunctionBucket( buffer.str(), this ) );
    (*coeff).coeff_fill(i);                                          // fill the coefficient (this won't do much for coefficient
                                                                     // functions)
    register_coeff(coeff, (*coeff).name());                          // register this coefficient in the system

  }

}

//*******************************************************************|************************************************************//
// fill in the data about each solver of this system 
//*******************************************************************|************************************************************//
void SpudSystemBucket::solvers_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  
  buffer.str("");  buffer << optionpath() << "/nonlinear_solver";
  int nsolvers = Spud::option_count(buffer.str());                   // find out how many nonlinear solvers there are
  for (uint i = 0; i < nsolvers; i++)                                // loop over them
  {
    buffer.str(""); buffer << optionpath() << "/nonlinear_solver[" 
                                                        << i << "]";

                                                                     // initialize a new solver bucket assuming this system is its
                                                                     // parent
    SpudSolverBucket_ptr solver( new SpudSolverBucket( buffer.str(), this ) );
    (*solver).fill();                                                // fill in the data about this solver bucket
    register_solver(solver, (*solver).name());                       // register the solver bucket in the system
  }
}

