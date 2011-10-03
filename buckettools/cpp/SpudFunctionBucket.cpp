
#include "BoostTypes.h"
#include "SpudFunctionBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemsWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudFunctionBucket::SpudFunctionBucket(const std::string &optionpath, 
                                            SystemBucket* system) : 
                                            optionpath_(optionpath), 
                                            FunctionBucket(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudFunctionBucket::~SpudFunctionBucket()
{
  empty_();                                                  // empty the data in the derived class
}

//*******************************************************************|************************************************************//
// fill the function bucket data structures assuming the buckettools schema and that this is a field
//*******************************************************************|************************************************************//
void SpudFunctionBucket::field_fill(const uint &index)
{
  std::stringstream buffer;                                          // optionpath buffer

  base_fill_(index);                                                 // fill the base data about the function bucket
                                                                     // (could be called from constructor as it's common regardless
                                                                     // of functionbucket type)

                                                                     // check we think this is a field (necessary but not
                                                                     // sufficient)
  assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field::"+name()));

  initialize_field_();                                               // initialize as a field

  functionals_fill_();                                               // fill in data about the functionals (doesn't attach
                                                                     // coefficients)

  if(include_in_visualization())
  {
                                                                     // initialize a pvd file (field specific so could be moved)
    pvdfile_.reset( new dolfin::File((*(*system_).bucket()).output_basename()+"_"+(*system_).name()+"_"+name()+".pvd", "compressed") );
  }

}

//*******************************************************************|************************************************************//
// fill the function bucket data structures assuming the buckettools schema and that this is a coefficient
//*******************************************************************|************************************************************//
void SpudFunctionBucket::coeff_fill(const uint &index)
{
  std::stringstream buffer;                                          // optionpath buffer

  base_fill_(index);                                                 // fill the base data about the function bucket
                                                                     // (called from constructor instead?)

  if(system_)                                                        // can't do this if this coefficient lives outside a system
  { 
                                                                     // check we think this is a coefficient (not a foolproof test)
    assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/coefficient::"+name()));
  }

  initialize_expression_coeff_();                                    // initialize as a coefficient (doesn't do much for coefficient
                                                                     // functions)

  functionals_fill_();                                               // fill in data about functionals (doesn't attach coefficients)

}

//*******************************************************************|************************************************************//
// initialize a coefficient that is of type function - left until almost last as we need a functionspace for it
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_function_coeff()
{
  std::stringstream buffer;                                          // optionpath buffer

  if (type_=="Function")                                             // only do this for coefficient functions
  {
    if(!(*(*system_).bucket()).contains_coefficientspace(uflsymbol()))// check we have a coefficient space for this function's ufl
    {                                                                // symbol
      dolfin::log(dolfin::ERROR, "Coefficient %s declared as type::Function but its ufl_symbol was not found in any forms or functionals.", 
                                 name().c_str());
      dolfin::error("Unable to allocate coefficient as no matching functionspace found.");
    }

    FunctionSpace_ptr coefficientspace =                             // grab the functionspace for this coefficient from the bucket
          (*(*system_).bucket()).fetch_coefficientspace(uflsymbol());// data maps

    function_.reset( new dolfin::Function(*coefficientspace) );      // allocate the function on this functionspace
    oldfunction_.reset( new dolfin::Function(*coefficientspace) );   // allocate the old function on this functionspace
    iteratedfunction_ = function_;                                   // just point this at the function

    buffer.str(""); buffer << optionpath() << "/type/rank/value";    // initialize the coefficient
    int nvs = Spud::option_count(buffer.str());
    if (nvs > 1)
    {
      dolfin::error("Haven't thought about values over regions.");
    }
    else
    {
//      for (uint i = 0; i < nvs; i++)
//      {
        uint i = 0;
        buffer.str(""); buffer << optionpath() << "/type/rank/value["
                                                         << i << "]";
        Expression_ptr valueexp = initialize_expression(buffer.str(),// initialize an expression from the optionpath 
                                                      &size_, &shape_);

                                                                     // interpolate this expression onto the function
        (*boost::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*valueexp);
        (*boost::dynamic_pointer_cast< dolfin::Function >(oldfunction_)).interpolate(*valueexp);

//      }
    }


  }

}

//*******************************************************************|************************************************************//
// make a partial copy of the provided function bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SpudFunctionBucket::copy_diagnostics(FunctionBucket_ptr &function, SystemBucket_ptr &system) const
{

  if(!function)
  {
    function.reset( new SpudFunctionBucket(optionpath_, &(*system)) );
  }

  FunctionBucket::copy_diagnostics(function, system);

  (*boost::dynamic_pointer_cast< SpudFunctionBucket >(function)).functional_optionpaths_ = functional_optionpaths_;

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a functional form in the function bucket data maps (with an optionpath as well)
//*******************************************************************|************************************************************//
void SpudFunctionBucket::register_functional(Form_ptr functional, 
                                      const std::string &name, 
                                      const std::string &optionpath)
{
  Form_it f_it = functionals_.find(name);                            // check if name already exists
  if (f_it != functionals_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
            "Functional named \"%s\" already exists in function.", 
                                                    name.c_str());
  }
  else
  {
    functionals_[name]            = functional;                      // it not, insert the form into the maps
    functional_optionpaths_[name] = optionpath;                      // and its optionpath too
  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in visualization output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_visualization() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_visualization");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in diagnostic output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_statistics() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_statistics");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in steady state check and output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_steadystate() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_steady_state");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in detectors output
//*******************************************************************|************************************************************//
const bool SpudFunctionBucket::include_in_detectors() const
{
  return Spud::have_option(optionpath()+"/diagnostics/include_in_detectors");
}

//*******************************************************************|************************************************************//
// return a string describing the contents of this function bucket
//*******************************************************************|************************************************************//
const std::string SpudFunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << " (" 
                                << optionpath() << ")" << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the functionals in this function bucket
//*******************************************************************|************************************************************//
const std::string SpudFunctionBucket::functionals_str(const int &indent) 
                                                              const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = functional_optionpaths_.begin(); 
                      s_it != functional_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Functional " << (*s_it).first << " (" 
                              << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// fill the function bucket base data assuming the buckettools schema (common for fields and coefficients)
//*******************************************************************|************************************************************//
void SpudFunctionBucket::base_fill_(const uint &index)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  index_ = index;                                                    // the field or coefficient index (most relevant for fields)

  buffer.str(""); buffer << optionpath() << "/name";                 // field or coefficient name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/ufl_symbol";           // ufl symbol
  serr = Spud::get_option(buffer.str(), uflsymbol_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/name";            // field or coefficient type (function, expression, constant)
  serr = Spud::get_option(buffer.str(), type_);                      // (as string)
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/rank/name";       // field or coefficient rank (as string)
  serr = Spud::get_option(buffer.str(), rank_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath()                             // the size of the element (if a vector)
                                      << "/type/rank/element/size";  // NOTE: defaults to the bucket geometry dimension
  serr = Spud::get_option(buffer.str(), size_, (*(*system_).bucket()).dimension()); 
  spud_err(buffer.str(), serr);

  std::vector< int > default_shape(2, (*(*system_).bucket()).dimension());
  buffer.str(""); buffer << optionpath()                             // the shape of the element (if a tensor)
                                    << "/type/rank/element/shape";   // NOTE: defaults to the bucket geometry dimension square
  serr = Spud::get_option(buffer.str(), shape_, default_shape); 
  spud_err(buffer.str(), serr);

  if (type_=="Constant")                                             // if the type is constant then we won't have element
  {                                                                  // information with which to set the size and shape
                                                                     // so they will automatically be defaults... check if that's
                                                                     // right
    buffer.str(""); buffer << optionpath() 
                                    << "/type/rank/value/constant";
    if (rank_=="Vector")                                             // get spud to tell us the element shape... if it's a vector
    {                                                                // the size is just the first entry of this
      std::vector< int > constant_shape;
      serr = Spud::get_option_shape(buffer.str(), constant_shape); 
      spud_err(buffer.str(), serr);
      size_ = constant_shape[0];
    }
    if (rank_=="Tensor")                                             // get the spud element shape for a tensor
    {
      serr = Spud::get_option_shape(buffer.str(), shape_); 
      spud_err(buffer.str(), serr);
    }
  }

}

//*******************************************************************|************************************************************//
// initialize a field assuming the base data has already been filled and using a buckettools spud schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_field_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  assert(type_=="Function");                                         // fields may only be functions

  if (rank_=="Tensor")                                               // broken (untested?) for field tensors
  {
    dolfin::error("Tensor fields not hooked up yet, sorry.");
  }
  
  int nfields =                                                      // how may fields in the system (i.e. mixed functionspace?)
      Spud::option_count((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field");
  if (nfields == 1)                                                  // no subfunctions (this field is identical to the system)
  {                                                                  // just grab references to:
    functionspace_ = (*system_).functionspace();                     // the functionspace

    function_ = (*system_).function();                               // the function
    oldfunction_ = (*system_).oldfunction();                         // the old function
    iteratedfunction_ = (*system_).iteratedfunction();               // the iterated function

    if (include_in_steadystate())
    {
      changefunction_ = (*system_).changefunction();                 // the change in the function between timesteps
    }
  }
  else                                                               // yes, multiple fields in this system so we need to make
  {                                                                  // a subspace and subfunctions (dangerous, be careful!)
    functionspace_ = (*(*system_).functionspace())[index_];          // create a subspace

    function_.reset( &(*(*system_).function())[index_],              // create a subfunction but don't delete it ever as it shares
                                              dolfin::NoDeleter() ); // information with the system function
    oldfunction_.reset( &(*(*system_).oldfunction())[index_],        // same for the old function
                                              dolfin::NoDeleter() );
    iteratedfunction_.reset( &(*(*system_).iteratedfunction())[index_], 
                                              dolfin::NoDeleter() ); // and the iterated function

    if (include_in_steadystate())
    {
      changefunction_.reset( &(*(*system_).changefunction())[index_],
                                              dolfin::NoDeleter() ); // and the change in the function between timesteps
    }
  }

  buffer.str(""); buffer << (*system_).name() << "::" << name();     // rename the function as SystemName::FieldName
  (*function_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << (*system_).name() << "::Old" << name();  // rename the old function as SystemName::OldFieldName
  (*oldfunction_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << (*system_).name() << "::Iterated"        // rename the iterated function as SystemName::IteratedFieldName
                                                         << name();  
  (*iteratedfunction_).rename(buffer.str(), buffer.str());

  if (include_in_steadystate())
  {
    buffer.str(""); buffer << (*system_).name() << "::TimestepChange"// rename the change function as SystemName::TimestepChangeFieldName
                                                           << name();  
    (*changefunction_).rename(buffer.str(), buffer.str());

    change_.reset( new double(0.0) );

    change_calculated_.reset( new bool(false) );

    buffer.str(""); buffer << optionpath() 
                      << "/diagnostics/include_in_steady_state/norm";
    serr = Spud::get_option(buffer.str(), change_normtype_);
    spud_err(buffer.str(), serr);
  }

  buffer.str(""); buffer << optionpath()                             // find out how many strong bcs there are (weak bcs handled in
                                << "/type/rank/boundary_condition";  // forms)
  int nbcs = Spud::option_count(buffer.str());
  if (nbcs > 0)                                                      // if we have any...
  {
    MeshFunction_uint_ptr edgeidmeshfunction =                       // get the edge id mesh function
      (*(*system_).mesh()).data().mesh_function("exterior_facet_domains");

    for (uint i = 0; i < nbcs; i++)                                  // loop over the bcs
    {
      buffer.str(""); buffer << optionpath() 
                    << "/type/rank/boundary_condition[" << i << "]";
      bc_fill_(buffer.str(), edgeidmeshfunction);                    // and fill in details about the bc
    }
  }
    
  buffer.str(""); buffer << optionpath()                             // find out how many initial conditions we have
                                  << "/type/rank/initial_condition";
  int nics = Spud::option_count(buffer.str());
  if (nics > 1)
  {
    dolfin::error("Haven't thought about ics over regions.");
  }
  else
  {
//    for (uint i = 0; i < nics; i++)
//    {
      uint i = 0;
      buffer.str(""); buffer << optionpath() 
                      << "/type/rank/initial_condition[" << i << "]";
      ic_fill_(buffer.str());                                        // fill in data about the field initial condition
//    }
  }
   
}

//*******************************************************************|************************************************************//
// fill the details of a bc assuming the buckettools schema and given a set of edge ids in a mesh function
//*******************************************************************|************************************************************//
void SpudFunctionBucket::bc_fill_(const std::string &optionpath,
                                  const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  
  std::string bcname;                                                // bc name
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), bcname); 
  spud_err(buffer.str(), serr);
  
  std::vector<int> bcids;                                            // the boundary ids this bc applies to
  buffer.str(""); buffer << optionpath << "/boundary_ids";
  serr = Spud::get_option(buffer.str(), bcids); 
  spud_err(buffer.str(), serr);
  
  buffer.str(""); buffer << optionpath << "/sub_components";         // how many subcomponents exist?
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsubcomp; i++)                                // loop over the subcomponents
  {
    buffer.str(""); buffer << optionpath << "/sub_components[" 
                                                        << i << "]";
    bc_component_fill_(buffer.str(), bcname, bcids,                  // initialize this component
                                      edgeidmeshfunction);
  }
}

//*******************************************************************|************************************************************//
// fill the details of a component of a bc on some boundary ids  assuming the buckettools schema and given a set of edge ids in a 
// mesh function
//*******************************************************************|************************************************************//
void SpudFunctionBucket::bc_component_fill_(const std::string &optionpath,
                                            const std::string &bcname,
                                            const std::vector<int> &bcids,
                                            const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  std::stringstream namebuffer;                                      // a buffer to assemble this bc's component name

  buffer.str(""); buffer << optionpath << "/components";             // do we have any components (i.e. could be scalar)?
  if (Spud::have_option(buffer.str()))
  {                                                                  // FIXME: tensorial components assumed broken!
    std::vector<int> subcompids;
    serr = Spud::get_option(buffer.str(), subcompids);               // get a list of the subcomponents 
    spud_err(buffer.str(), serr);
    
    for (std::vector<int>::const_iterator subcompid =                // loop over those subcomponents
                                          subcompids.begin(); 
                      subcompid < subcompids.end(); subcompid++)
    {

       FunctionSpace_ptr subfunctionspace =                          // grab the subspace from the field functionspace
                                    (*functionspace())[*subcompid];  // happily dolfin caches this for us
       
       buffer.str(""); buffer << optionpath << "/type::Dirichlet";   // initialize an expression describing the function
       Expression_ptr bcexp =                                        // on this boundary condition, assuming it is of type
                  initialize_expression(buffer.str(), &size_, &shape_);// Dirichlet
       
       namebuffer.str(""); namebuffer << bcname << "::" 
                                                       << *subcompid;// assemble a name for the bc (including the subcomponent id)
       register_bcexpression(bcexp, namebuffer.str());               // register the expression in the function bucket maps
       
       for (std::vector<int>::const_iterator bcid = bcids.begin();   // loop over the boundary ids
                                          bcid < bcids.end(); bcid++)
       {                                                             // create a new bc on each boundary id for this subcomponent
         BoundaryCondition_ptr bc(new dolfin::DirichletBC(*subfunctionspace, *bcexp, *edgeidmeshfunction, *bcid));
         namebuffer.str(""); namebuffer << bcname << "::" 
                                      << *subcompid << "::" << *bcid;// assemble a name incorporating the boundary id
         register_bc(bc, namebuffer.str());                          // register the bc in the function bucket
       }
       
    }
  }
  else                                                               // no components (scalar or using all components)
  {
    buffer.str(""); buffer << optionpath << "/type::Dirichlet";      // initialize an expression describing the field on the
    Expression_ptr bcexp =                                           // boundary
                  initialize_expression(buffer.str(), &size_, &shape_);
    
    register_bcexpression(bcexp, bcname);                            // register the expression in the function bucket
    
    for (std::vector<int>::const_iterator bcid = bcids.begin();      // loop over the boundary ids
                                          bcid < bcids.end(); bcid++)
    {                                                                // create a bc on each boundary id for all components
      BoundaryCondition_ptr bc(new dolfin::DirichletBC(*functionspace(), *bcexp, *edgeidmeshfunction, *bcid));
      namebuffer.str(""); namebuffer << bcname << "::" << *bcid;     // assemble a name for this bc incorporating the boundary id
      register_bc(bc, namebuffer.str());                             // register the bc in the function bucket
    }
  }
}

//*******************************************************************|************************************************************//
// fill the details of an initial condition assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::ic_fill_(const std::string &optionpath)
{
  icexpression_ = initialize_expression(optionpath, &size_, &shape_);  // intialize the expression (FIXME: still needs regions support!)
}

//*******************************************************************|************************************************************//
// initialize a coefficient that is of type expression or constant
//*******************************************************************|************************************************************//
void SpudFunctionBucket::initialize_expression_coeff_()
{
  std::stringstream buffer;                                          // optionpath buffer

  if (type_=="Expression" || type_=="Constant")                      // only proceed for expressions and constants
  {

    buffer.str(""); buffer << optionpath() << "/type/rank/value";    // find out how many value regions we have
    int nvs = Spud::option_count(buffer.str());
    if (nvs > 1)
    {
      dolfin::error("Haven't thought about values over regions.");
    }
    else
    {
//      for (uint i = 0; i < nvs; i++)
//      {
        uint i = 0;
        buffer.str(""); buffer << optionpath() << "/type/rank/value[" 
                                                        << i << "]";

        function_ = 
                initialize_expression(buffer.str(), &size_, &shape_);// initialize the function from the optionpath
                                                                     // (this will do almost nothing for functionals)
        oldfunction_ = 
                initialize_expression(buffer.str(), &size_, &shape_);// initialize the function from the optionpath
                                                                     // (this will do almost nothing for functionals)

        iteratedfunction_ = function_;                               // for now just point this at the function_

        buffer.str(""); buffer << optionpath() << "/type/rank/value[" 
                                              << i << "]/functional";
        if(Spud::have_option(buffer.str()))
        {
          constantfunctional_fill_();
        }

//      }
    }
  }

}

//*******************************************************************|************************************************************//
// fill the details of a constant functional assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::constantfunctional_fill_()
{
  constantfunctional_ = ufc_fetch_functional((*system_).name(),      // get a pointer to the functional form from the ufc
                        name(), (*system_).mesh());

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // functional because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
  uint ncoeff = (*constantfunctional_).num_coefficients();           // how many coefficients does this functional require?
  for (uint i = 0; i < ncoeff; i++)
  {
    std::string uflsymbol = 
                          (*constantfunctional_).coefficient_name(i);// what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
    if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))    // a base ufl symbol was only inserted into the parent bucket's
    {                                                                // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
      std::string baseuflsymbol =                                    // what is the base ufl symbol?
            (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol);   // have we already registered a functionspace for this base ufl
                                                                     // symbol?
      if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
      {                                                              // no...
        FunctionSpace_ptr coefficientspace;
        coefficientspace = 
                         ufc_fetch_coefficientspace_from_functional( // take a pointer to the functionspace from the ufc
                                      (*system_).name(), name(), 
                                      baseuflsymbol, 
                                      (*system_).mesh());
        (*(*system_).bucket()).register_coefficientspace(            // and register it in the parent bucket's map
                                      coefficientspace, 
                                      baseuflsymbol);
      }
    }

  }
}

//*******************************************************************|************************************************************//
// fill the details of any functionals associated with this field or coefficient assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudFunctionBucket::functionals_fill_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
   
  buffer.str(""); buffer << optionpath() 
                     << "/diagnostics/include_in_statistics/functional";
  int nfuncs = Spud::option_count(buffer.str());                     // find out how many functionals there are
  for (uint i = 0; i < nfuncs; i++)
  {
    buffer.str(""); buffer << optionpath()  
        << "/diagnostics/include_in_statistics/functional[" << i << "]";
    std::string functionaloptionpath = buffer.str();

    std::string functionalname;                                      // the functional name
    buffer.str(""); buffer << functionaloptionpath << "/name";
    serr = Spud::get_option(buffer.str(), functionalname); 
    spud_err(buffer.str(), serr);

    Form_ptr functional = ufc_fetch_functional((*system_).name(),    // get a pointer to the functional form from the ufc
                          name(), functionalname, (*system_).mesh());
    register_functional(functional, functionalname,                  // register it in the function bucket
                                               functionaloptionpath);

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // functional because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
    uint ncoeff = (*functional).num_coefficients();                  // how many coefficients does this functional require?
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*functional).coefficient_name(i);     // what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
      if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))  // a base ufl symbol was only inserted into the parent bucket's
      {                                                              // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
        std::string baseuflsymbol =                                  // what is the base ufl symbol?
              (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol); // have we already registered a functionspace for this base ufl
                                                                     // symbol?
        if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
        {                                                            // no...
          FunctionSpace_ptr coefficientspace;
          coefficientspace = 
                        ufc_fetch_coefficientspace_from_functional(  // take a pointer to the functionspace from the ufc
                                        (*system_).name(), name(), 
                                        functionalname, baseuflsymbol, 
                                        (*system_).mesh());
          (*(*system_).bucket()).register_coefficientspace(          // and register it in the parent bucket's map
                                        coefficientspace, 
                                        baseuflsymbol);
        }
      }

    }
  }

}

//*******************************************************************|************************************************************//
// empty the data structures in the spudfunctionbucket
//*******************************************************************|************************************************************//
void SpudFunctionBucket::empty_()
{
  functional_optionpaths_.clear();
  FunctionBucket::empty_();
}

