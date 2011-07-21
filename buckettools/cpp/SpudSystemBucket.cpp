
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

// Specific constructor
SpudSystemBucket::SpudSystemBucket(std::string optionpath, Bucket* bucket) : optionpath_(optionpath), SystemBucket(bucket)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudSystemBucket::~SpudSystemBucket()
{
  // Do nothing
}

// Fill the system using spud and assuming a buckettools schema structure
void SpudSystemBucket::fill()
{
  // fill in the base data (could be called from the constructor)
  base_fill_();

  // register the functionspace and (mixed?) functions associated with this system
  systemfunction_fill_();

  // Having just registered the system uflsymbols, quickly loop through
  // the rest of the fields recording their symbol and a null pointer
  uflsymbols_fill_();

  // initialize the fields
  fields_fill_(); 
 
  // initialize the solvers
  solvers_fill_();

  // initialize the coefficients
  // (we do this after the solvers so that we can get the coefficientspaces from the solver forms)
  coeffs_fill_(); 

  // attach the coefficients to the functionals and forms
  // this will probably have to move to later (after aliased coefficients are populated)
  attach_all_coeffs_();
}

void SpudSystemBucket::base_fill_()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  Spud::OptionError serr;

  // Get the system name
  buffer.str(""); buffer << optionpath() << "/name";
  serr = Spud::get_option(buffer.str(), name_); spud_err(buffer.str(), serr);

  // Get the coeff ufl symbol (necessary to register the field name)
  buffer.str(""); buffer << optionpath() << "/ufl_symbol";
  serr = Spud::get_option(buffer.str(), uflsymbol_); spud_err(buffer.str(), serr);

  // Get the name of the mesh this system is defined on
  std::string meshname;
  buffer.str(""); buffer << optionpath() << "/mesh/name";
  serr = Spud::get_option(buffer.str(), meshname); spud_err(buffer.str(), serr);
  // and then extract the mesh from the bucket we're filling
  mesh_ = (*bucket_).fetch_mesh(meshname);

}

void SpudSystemBucket::systemfunction_fill_()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  Spud::OptionError serr;

  // Here's where the automatically generated magic happens... this is fetching
  // the functionspace from ufl
  functionspace_ = ufc_fetch_functionspace(name(), mesh());

  // Declare a series of functions on this functionspace:
  // Current
  function_.reset( new dolfin::Function(*functionspace_) );
  buffer.str(""); buffer << name() << "::Function";
  (*function_).rename( buffer.str(), buffer.str() );

  // Old
  oldfunction_.reset( new dolfin::Function(*functionspace_) );
  buffer.str(""); buffer << name() << "::OldFunction";
  (*oldfunction_).rename( buffer.str(), buffer.str() );

  // Iterated
  iteratedfunction_.reset( new dolfin::Function(*functionspace_) );
  buffer.str(""); buffer << name() << "::IteratedFunction";
  (*iteratedfunction_).rename( buffer.str(), buffer.str() );

  register_uflsymbol(function_, uflsymbol_);
  register_uflsymbol(oldfunction_, uflsymbol_+"_n");
  register_uflsymbol(iteratedfunction_, uflsymbol_+"_i");
  // we intentionally do not register the ufl with the function name
  // as we only do that for functions not systems (in case they share
  // a name)

}

void SpudSystemBucket::uflsymbols_fill_()
{
  std::stringstream buffer;
  Spud::OptionError serr;
  std::string uflsymbol, functionname;

  buffer.str("");  buffer << optionpath() << "/field";
  int nfields = Spud::option_count(buffer.str());
  // Loop over the fields (which are subfunctions of this functionspace)
  // and register their uflsymbols in the system
  for (uint i = 0; i < nfields; i++)
  {
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]/name";
    serr = Spud::get_option(buffer.str(), functionname); spud_err(buffer.str(), serr);
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]/ufl_symbol";
    serr = Spud::get_option(buffer.str(), uflsymbol); spud_err(buffer.str(), serr);
    register_uflname(functionname, uflsymbol);
    register_uflname(functionname, uflsymbol+"_i");
    register_uflname(functionname, uflsymbol+"_n");
  }

  buffer.str("");  buffer << optionpath() << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());
  // Loop over the coefficients and register their ufsymbols
  for (uint i = 0; i < ncoeffs; i++)
  {
    buffer.str(""); buffer << optionpath() << "/coefficient[" << i << "]/name";
    serr = Spud::get_option(buffer.str(), functionname); spud_err(buffer.str(), serr);
    buffer.str(""); buffer << optionpath() << "/coefficient[" << i << "]/ufl_symbol";
    serr = Spud::get_option(buffer.str(), uflsymbol); spud_err(buffer.str(), serr);
    register_uflname(functionname, uflsymbol);
    register_uflname(functionname, uflsymbol+"_i");
    register_uflname(functionname, uflsymbol+"_n");
  }
}

// Fill out the information regarding the each subfunction (or field)
void SpudSystemBucket::fields_fill_()
{
  std::stringstream buffer;

  // A counter for the components in this system (allows the ic to be generalized)
  uint component = 0;
  std::map< uint, Expression_ptr > icexpressions;

  buffer.str("");  buffer << optionpath() << "/field";
  int nfields = Spud::option_count(buffer.str());
  // Loop over the fields (which are subfunctions of this functionspace)
  // and register them in the system
  for (uint i = 0; i < nfields; i++)
  {
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]";

    SpudFunctionBucket_ptr field( new SpudFunctionBucket( buffer.str(), this ) );
    (*field).field_fill(i);
    register_field(field, (*field).name());

    // set the function associated with the uflsymbol back in the system
    register_uflsymbol((*field).function(), (*field).uflsymbol());
    register_uflsymbol((*field).oldfunction(), (*field).uflsymbol()+"_n");
    register_uflsymbol((*field).iteratedfunction(), (*field).uflsymbol()+"_i");

    uint_Expression_it e_it = icexpressions.find(component);
    if (e_it != icexpressions.end())
    {
      dolfin::error("ICExpression with component number %d already exists in icexpressions map.", component);
    }
    else
    {
      icexpressions[component] = (*field).icexpression();
    }

    component += (*(*field).icexpression()).value_size();
  }

  // While filling the fields we should have set up a map from
  // components to initial condition expressions... use this
  // now to initialize the whole system function to that
  // ic
  apply_ic_(component, icexpressions);

}

void SpudSystemBucket::apply_ic_(const uint &component, const std::map< uint, Expression_ptr > &icexpressions)
{
  Expression_ptr ic;
  if (component==1)
  {
    ic.reset( new InitialConditionExpression(icexpressions) );
  }
  else
  {
    ic.reset( new InitialConditionExpression(component, icexpressions));
  }
  (*oldfunction_).interpolate(*ic);
  *iteratedfunction_ = *oldfunction_;
  *function_ = *oldfunction_;
}

//void SpudSystemBucket::apply_bc_()
//{
//  for(std::map< std::string, DirichletBC_ptr >::iterator
//            bc = dirichletbcs_begin(); 
//            bc != dirichletbcs_end(); bc++)
//  {
//    (*((*bc).second)).apply((*sysfunc).vector());
//  }
//  
//}

// Fill out the information regarding each coefficient
void SpudSystemBucket::coeffs_fill_()
{
  std::stringstream buffer;
  
  buffer.str("");  buffer << optionpath() << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());
  // Loop over the coefficients and register those that are expressions or constants
  // Can't do much with the functions because they depend on solvers
  // Aliased coefficients depend on all the other systems being populated first
  for (uint i = 0; i < ncoeffs; i++)
  {
    buffer.str(""); buffer << optionpath() << "/coefficient[" << i << "]";

    SpudFunctionBucket_ptr coeff( new SpudFunctionBucket( buffer.str(), this ) );
    (*coeff).coeff_fill(i);
    register_coeff(coeff, (*coeff).name());

  }

  // std::cout << coefficientspaces_str() << std::flush;

  // Go for a second loop over the coefficients as we didn't initialize functions last time
  for (FunctionBucket_it f_it = coeffs_.begin(); f_it != coeffs_.end(); f_it++)
  {
    if ((*(*f_it).second).type()=="Function")
    {
      (*(boost::dynamic_pointer_cast< SpudFunctionBucket >((*f_it).second))).initialize_function_coeff();
    }

    // reset the function associated with the uflsymbol back in the system
    register_uflsymbol((*(*f_it).second).function(), (*(*f_it).second).uflsymbol());
    register_uflsymbol((*(*f_it).second).oldfunction(), (*(*f_it).second).uflsymbol()+"_n");
    register_uflsymbol((*(*f_it).second).iteratedfunction(), (*(*f_it).second).uflsymbol()+"_i");
  }

}

void SpudSystemBucket::solvers_fill_()
{
  std::stringstream buffer;
  
  buffer.str("");  buffer << optionpath() << "/nonlinear_solver";
  int nsolvers = Spud::option_count(buffer.str());
  // Loop over the nonlinear_solvers and register them
  for (uint i = 0; i < nsolvers; i++)
  {
    buffer.str(""); buffer << optionpath() << "/nonlinear_solver[" << i << "]";

    SpudSolverBucket_ptr solver( new SpudSolverBucket( buffer.str(), this ) );
    (*solver).fill();
    register_solver(solver, (*solver).name());
  }
}

// Return a string describing the contents of the system
std::string SpudSystemBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SystemBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << uflsymbols_str(indent);
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << solvers_str(indent);
  return s.str();
}

