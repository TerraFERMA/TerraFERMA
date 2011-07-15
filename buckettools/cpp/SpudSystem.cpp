
#include "BoostTypes.h"
#include "InitialConditionExpression.h"
#include "SpudSystem.h"
#include "System.h"
#include "SystemsWrapper.h"
#include "SpudBase.h"
#include "SpudFunctionBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

// Specific constructor
SpudSystem::SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh, Bucket* bucket) : optionpath_(optionpath), System(name, mesh, bucket)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudSystem::~SpudSystem()
{
  // Do nothing
}

// Fill the system using spud and assuming a buckettools schema structure
void SpudSystem::fill()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;

  // Here's where the automatically generated magic happens... this is fetching
  // the functionspace from ufl
  functionspace_ = ufc_fetch_functionspace(name(), mesh_);

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

  // A counter for the components in this system (allows the ic to be generalized)
  uint component = 0;
  buffer.str("");  buffer << optionpath() << "/field";
  int nfields = Spud::option_count(buffer.str());
  // Loop over the fields (which are subfunctions of this functionspace)
  // and register them in the system
  for (uint i = 0; i < nfields; i++)
  {
    buffer.str(""); buffer << optionpath() << "/field[" << i << "]";
    fields_fill_(buffer.str(), i, nfields, component); 
  }

  // While filling the fields we should have set up a map from
  // components to initial condition expressions... use this
  // now to initialize the whole system function to that
  // ic
  Expression_ptr ic;
  if (component==1)
  {
    ic.reset( new InitialConditionExpression(icexpressions_) );
  }
  else
  {
    ic.reset( new InitialConditionExpression(component, icexpressions_));
  }
  (*oldfunction_).interpolate(*ic);
  *iteratedfunction_ = *oldfunction_;
  *function_ = *oldfunction_;
//  for(std::map< std::string, DirichletBC_ptr >::iterator
//            bc = dirichletbcs_begin(); 
//            bc != dirichletbcs_end(); bc++)
//  {
//    (*((*bc).second)).apply((*sysfunc).vector());
//  }
//  

  buffer.str("");  buffer << optionpath() << "/coefficient";
  int ncoeffs = Spud::option_count(buffer.str());
  // Loop over the coefficients and register those that are expressions or constants
  // Can't do anything with the functions because they depend on solvers
  // Aliased coefficients depend on all the other systems being populated first
  for (uint i = 0; i < ncoeffs; i++)
  {
    buffer.str(""); buffer << optionpath() << "/coefficient[" << i << "]";
    coeffs_fill_(buffer.str()); 
  }

}

// Fill out the information regarding the each subfunction (or field)
void SpudSystem::fields_fill_(const std::string &optionpath, 
                              const uint &field_i, 
                              const uint &nfields, 
                              uint &component)
{
  // Fields are strange because lots of their data actually belongs
  // to the system so we do a lot of the population of the data structures
  // at a system level rather than at the function class level.
  // The function class (which will also get populated from here)
  // mostly exists to map the functionals for the stat output but these
  // don't get completed until later since they can depend on anything
  // in the system.

  // A buffer to put option paths in
  std::stringstream buffer;

  // Get the field type
  std::string type;
  buffer.str(""); buffer << optionpath << "/type/name";
  Spud::get_option(buffer.str(), type);
  // we only know how to deal with Functions at the moment but this may change
  assert(type=="Function");

  // What is the field size (if it's a vector)
  // Would it be possible to get this from the subfunctionspace below?
  int size;
  buffer.str(""); buffer << optionpath << "/type/rank/element/size";
  Spud::get_option(buffer.str(), size, (*bucket_).dimension());

  // What is the field shape (if it's a tensor)
  // Would it be possible to get this from the subfunctionspace below?
  std::vector< int > shape;
  std::vector< int > default_shape(2, (*bucket_).dimension());
  buffer.str(""); buffer << optionpath << "/type/rank/element/shape";
  Spud::get_option(buffer.str(), shape, default_shape);

  // Most of what is below is broken for tensors so let's just die here with an error message
  buffer.str(""); buffer << optionpath << "/type/rank::Tensor";
  if (Spud::have_option(buffer.str()))
  {
    dolfin::error("Tensor fields not hooked up yet, sorry.");
  }
  
  // Get the coeff ufl symbol (necessary to register the field name)
  std::string uflsymbol;
  buffer.str(""); buffer << optionpath << "/ufl_symbol";
  Spud::get_option(buffer.str(), uflsymbol);

  // What is this field called?
  std::string fieldname;
  buffer.str(""); buffer << optionpath << "/name";
  Spud::get_option(buffer.str(), fieldname);

  // Register this fielname with the unique (on a system level) uflsymbol
  register_uflname(fieldname, uflsymbol);

  // Is this a mixed functionspace or not?
  FunctionSpace_ptr subfunctionspace;
  Function_ptr function, oldfunction, iteratedfunction;
  if (nfields == 1)
  {
    // no... the subfunctionspace for this field is identical to the system's
    // luckily for us these are just pointers so grab a reference to it
    subfunctionspace = functionspace_;

    // not sure quite what this will do (in the nfields==1 case) but let's try to register the field
    function.reset( new dolfin::Function( (*function_) ) );
    oldfunction.reset( new dolfin::Function( (*oldfunction_) ) );
    iteratedfunction.reset( new dolfin::Function( (*iteratedfunction_) ) );
  }
  else
  {
    // yes... use DOLFIN to extract that subspace so we can declare things on it (ics, bcs etc.)
    subfunctionspace.reset( new dolfin::SubSpace(*functionspace_, field_i) );

    // not sure quite what this will do (in the nfields==1 case) but let's try to register the field
    function.reset( new dolfin::Function( (*function_)[field_i] ) );
    oldfunction.reset( new dolfin::Function( (*oldfunction_)[field_i] ) );
    iteratedfunction.reset( new dolfin::Function( (*iteratedfunction_)[field_i] ) );
  }
  // register a pointer to the subfunctionspace in the system so it remains in memory for other
  // objects that take a reference
  register_subfunctionspace(subfunctionspace, fieldname);
  // register a pointer to the field as well (first give it a sensible name and label)
  buffer.str(""); buffer << name() << "::" << fieldname;
  (*function).rename(buffer.str(), buffer.str());
  SpudFunctionBucket_ptr field( new SpudFunctionBucket( fieldname, optionpath, function, oldfunction, iteratedfunction, this ) );
  (*field).fill();
  register_field(field, fieldname, optionpath);

  buffer.str(""); buffer << optionpath << "/type/rank/boundary_condition";
  int nbcs = Spud::option_count(buffer.str());
  if (nbcs > 0)
  {
    // get the edge id information to set the bcs
    MeshFunction_uint_ptr edgeidmeshfunction = (*mesh_).data().mesh_function("EdgeIDs");

    for (uint i = 0; i < nbcs; i++)
    {
      buffer.str(""); buffer << optionpath << "/type/rank/boundary_condition[" << i << "]";
      bc_fill_(buffer.str(), fieldname, size, shape, subfunctionspace, edgeidmeshfunction);
    }
  }
    
  buffer.str(""); buffer << optionpath << "/type/rank/initial_condition";
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
      buffer.str(""); buffer << optionpath << "/type/rank/initial_condition[" << i << "]";
      ic_fill_(buffer.str(), size, shape, component);
//    }
  }
   
}

// Fill out the information regarding the each subfunction (or field)
void SpudSystem::coeffs_fill_(const std::string &optionpath)
{
  // A buffer to put option paths in
  std::stringstream buffer;

  // Get the coeff type
  std::string type;
  buffer.str(""); buffer << optionpath << "/type/name";
  Spud::get_option(buffer.str(), type);

  if (type=="Function" || type=="Aliased")
  {
    // These don't get handled here
  }
  else if (type=="Expression" || type=="Constant")
  {

    // Get the coeff ufl symbol (necessary for all sorts of things below)
    std::string uflsymbol;
    buffer.str(""); buffer << optionpath << "/ufl_symbol";
    Spud::get_option(buffer.str(), uflsymbol);

    // What is this coeff called?
    std::string coeffname;
    buffer.str(""); buffer << optionpath << "/name";
    Spud::get_option(buffer.str(), coeffname);

    register_uflname(coeffname, uflsymbol);

    // Most of what is below is broken for tensors so let's just die here with an error message
    std::string rank;
    buffer.str(""); buffer << optionpath << "/type/rank/name";
    Spud::get_option(buffer.str(), rank);
    if (rank=="Tensor")
    {
      dolfin::error("Tensor coefficients not hooked up yet, sorry.");
    }
  
    // What is the field size (if it's a vector)
    int size;
    buffer.str(""); buffer << optionpath << "/type/rank/element/size";
    Spud::get_option(buffer.str(), size, (*bucket_).dimension());

    // What is the field shape (if it's a tensor)
    std::vector< int > shape;
    std::vector< int > default_shape(2, (*bucket_).dimension());
    buffer.str(""); buffer << optionpath << "/type/rank/element/shape";
    Spud::get_option(buffer.str(), shape, default_shape);
    if (type=="Constant")
    {
      // These will have just been set to defaults in this case, which may not be right
      buffer.str(""); buffer << optionpath << "/type/value/constant";
      if (rank=="Vector")
      {
        std::vector< int > constant_shape;
        Spud::get_option_shape(buffer.str(), constant_shape);
        size = constant_shape[0];
      }
      if (rank=="Tensor")
      {
        Spud::get_option_shape(buffer.str(), shape);
      }
    }

    buffer.str(""); buffer << optionpath << "/type/rank/value";
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
        buffer.str(""); buffer << optionpath << "/type/rank/value[" << i << "]";
        coeff_exp_fill_(buffer.str(), size, shape, coeffname);
//      }
    }
  }
  else
  {
    dolfin::error("Unknown coeff type in coeffs_fill_.");
  }

}

void SpudSystem::bc_fill_(const std::string &optionpath,
                          const std::string &fieldname,
                          const int &size,
                          const std::vector<int> &shape,
                          const FunctionSpace_ptr &subfunctionspace, 
                          const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;
  
  std::string bcname;
  buffer.str(""); buffer << optionpath << "/name";
  Spud::get_option(buffer.str(), bcname);
  
  std::vector<int> bcids;
  buffer.str(""); buffer << optionpath << "/boundary_ids";
  Spud::get_option(buffer.str(), bcids);
  
  buffer.str(""); buffer << optionpath << "/sub_components";
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsubcomp; i++)
  {
    buffer.str(""); buffer << optionpath << "/sub_components[" << i << "]";
    bc_component_fill_(buffer.str(), fieldname, bcname, size, shape, bcids, subfunctionspace, edgeidmeshfunction);
  }
}

void SpudSystem::bc_component_fill_(const std::string &optionpath,
                        const std::string &fieldname,
                        const std::string &bcname,
                        const int &size,
                        const std::vector<int> &shape,
                        const std::vector<int> &bcids,
                        const FunctionSpace_ptr &subfunctionspace,
                        const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;
  std::stringstream namebuffer;

  buffer.str(""); buffer << optionpath << "/components";
  if (Spud::have_option(buffer.str()))
  {
    // FIXME: tensor support needs to go in here, so another switch between a list and a python function
    //        to return the components!
    std::vector<int> subcompids;
    Spud::get_option(buffer.str(), subcompids);
    
    for (std::vector<int>::const_iterator subcompid = subcompids.begin(); subcompid < subcompids.end(); subcompid++)
    {

       FunctionSpace_ptr subsubfunctionspace;
       namebuffer.str(""); buffer << fieldname << "::" << *subcompid;
       if (contains_subfunctionspace(namebuffer.str()))
       {
         subsubfunctionspace = fetch_subfunctionspace(namebuffer.str());
       }
       else
       {
         subsubfunctionspace.reset( new dolfin::SubSpace(*subfunctionspace, *subcompid) );
         register_subfunctionspace(subsubfunctionspace, buffer.str());
       }
       
       buffer.str(""); buffer << optionpath << "/type::Dirichlet";
       Expression_ptr bcexp = initialize_expression(buffer.str(), size, shape);
       
       namebuffer.str(""); namebuffer << fieldname << "::" << *subcompid << "::" << bcname;
       register_bcexpression(bcexp, namebuffer.str());
       
       for (std::vector<int>::const_iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
       {
         DirichletBC_ptr bc(new dolfin::DirichletBC(*subsubfunctionspace, *bcexp, *edgeidmeshfunction, *bcid));
         namebuffer.str(""); namebuffer << fieldname << "::" << *subcompid << "::" << bcname << "::" << *bcid;
         register_dirichletbc(bc, namebuffer.str());
       }
       
    }
  }
  else
  {
    buffer.str(""); buffer << optionpath << "/type::Dirichlet";
    Expression_ptr bcexp = initialize_expression(buffer.str(), size, shape);
    
    namebuffer.str(""); namebuffer << fieldname << "::" << bcname;
    register_bcexpression(bcexp, namebuffer.str());
    
    for (std::vector<int>::const_iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
    {
      DirichletBC_ptr bc(new dolfin::DirichletBC(*subfunctionspace, *bcexp, *edgeidmeshfunction, *bcid));
      namebuffer.str(""); namebuffer << fieldname << "::" << bcname << "::" << *bcid;
      register_dirichletbc(bc, namebuffer.str());
    }
  }
}

void SpudSystem::coeff_exp_fill_(const std::string &optionpath,
                                 const int &size,
                                 const std::vector<int> &shape,
                                 const std::string &coeffname)
{
  // This will get more complicated when we support looping over regions
  // Loop will probably have to move in here!
  Expression_ptr exp = initialize_expression(optionpath, size, shape);

  SpudFunctionBucket_ptr coeff( new SpudFunctionBucket( coeffname, optionpath, exp, this ) );
  (*coeff).fill();

  register_coeff(coeff, coeffname);
}

void SpudSystem::ic_fill_(const std::string &optionpath,
                          const int &size,
                          const std::vector<int> &shape,
                          uint &component)
{
  // This will get more complicated when we support looping over regions
  // Loop will probably have to move in here!
  Expression_ptr icexp = initialize_expression(optionpath, size, shape);

  register_icexpression(icexp, component);

  component += (*icexp).value_size();
}

// Register a dolfin function as a field in the system
void SpudSystem::register_field(FunctionBucket_ptr field, std::string name, std::string optionpath)
{
  // First check if a field with this name already exists
  FunctionBucket_it f_it = fields_.find(name);
  if (f_it != fields_.end())
  {
    // if it does, issue an error
    dolfin::error("Field named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    fields_[name]            = field;
    field_optionpaths_[name] = optionpath;
  }
}

// Return a string describing the contents of the system
std::string SpudSystem::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "System " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << fields_str(indent);
  s << coeffs_str(indent);
  s << bcexpressions_str(indent);
  s << icexpressions_str(indent);
  return s.str();
}

