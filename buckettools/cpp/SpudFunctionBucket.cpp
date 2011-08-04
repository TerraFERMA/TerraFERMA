
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

// Specific constructor
SpudFunctionBucket::SpudFunctionBucket(std::string optionpath, SystemBucket* system) : optionpath_(optionpath), FunctionBucket(system)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudFunctionBucket::~SpudFunctionBucket()
{
  // Do nothing
}

// Fill the function using spud and assuming a buckettools schema structure
void SpudFunctionBucket::field_fill(const uint &index)
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;

  base_fill_(index);

  // we do this first in case this is a function coefficient that's only referenced in this functional
  // - hence we'll have the functionspace available later
  functionals_fill_();

  assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field::"+name()));

  // this is a field
  initialize_field_();

  pvdfile_.reset( new dolfin::File((*system_).name()+"_"+name()+".pvd", "compressed") );

}

// Fill the function using spud and assuming a buckettools schema structure
void SpudFunctionBucket::coeff_fill(const uint &index)
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;

  base_fill_(index);

  // we do this first in case this is a function coefficient that's only referenced in this functional
  // - hence we'll have the functionspace available later
  functionals_fill_();

  assert(Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/coefficient::"+name()));

  // this is a coefficient
  initialize_expression_coeff_();

}

void SpudFunctionBucket::base_fill_(const uint &index)
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;

  // the index, the name, the optionpath, the ufl_symbol and the type
  // are the pieces of information unique to each field

  // the index of this field
  index_ = index;

  // What is this field called?
  buffer.str(""); buffer << optionpath() << "/name";
  serr = Spud::get_option(buffer.str(), name_); spud_err(buffer.str(), serr);

  // Get the coeff ufl symbol (necessary to register the field name)
  buffer.str(""); buffer << optionpath() << "/ufl_symbol";
  serr = Spud::get_option(buffer.str(), uflsymbol_); spud_err(buffer.str(), serr);

  // Get the field type
  buffer.str(""); buffer << optionpath() << "/type/name";
  serr = Spud::get_option(buffer.str(), type_); spud_err(buffer.str(), serr);

  // Most of what is below is broken for tensors so let's just die here with an error message
  buffer.str(""); buffer << optionpath() << "/type/rank/name";
  serr = Spud::get_option(buffer.str(), rank_); spud_err(buffer.str(), serr);

  // What is the function size (if it's a vector)
  // Would it be possible to get this from the subfunctionspace below?
  buffer.str(""); buffer << optionpath() << "/type/rank/element/size";
  serr = Spud::get_option(buffer.str(), size_, (*(*system_).bucket()).dimension()); spud_err(buffer.str(), serr);

  // What is the function shape (if it's a tensor)
  // Would it be possible to get this from the subfunctionspace below?
  std::vector< int > default_shape(2, (*(*system_).bucket()).dimension());
  buffer.str(""); buffer << optionpath() << "/type/rank/element/shape";
  serr = Spud::get_option(buffer.str(), shape_, default_shape); spud_err(buffer.str(), serr);

  if (type_=="Constant")
  {
    // These will have just been set to defaults in this case, which may not be right
    // - only coefficients can end up in here so only consider their optionpaths
    buffer.str(""); buffer << optionpath() << "/type/rank/value/constant";
    if (rank_=="Vector")
    {
      std::vector< int > constant_shape;
      serr = Spud::get_option_shape(buffer.str(), constant_shape); spud_err(buffer.str(), serr);
      size_ = constant_shape[0];
    }
    if (rank_=="Tensor")
    {
      serr = Spud::get_option_shape(buffer.str(), shape_); spud_err(buffer.str(), serr);
    }
  }

}

void SpudFunctionBucket::initialize_field_()
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;

  // we only know how to deal with Functions at the moment but this may change
  assert(type_=="Function");

  // Most of what is below is broken for tensors so let's just die here with an error message
  if (rank_=="Tensor")
  {
    dolfin::error("Tensor fields not hooked up yet, sorry.");
  }
  
  // Is this a mixed functionspace or not (i.e. how many fields does the system have)?
  int nfields = Spud::option_count((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/field");
  if (nfields == 1)
  {
    // no... the subfunctionspace for this field is identical to the system's
    // luckily for us these are just pointers so grab a reference to it
    functionspace_ = (*system_).functionspace();

    // not sure quite what this will do (in the nfields==1 case) but let's try to register the field
    function_ = (*system_).function();
    oldfunction_ = (*system_).oldfunction();
    iteratedfunction_ = (*system_).iteratedfunction();
  }
  else
  {
    // yes... use DOLFIN to extract that subspace so we can declare things on it (ics, bcs etc.)
    functionspace_ = (*(*system_).functionspace())[index_];

    // register a shallow copy of the subfunction and make sure it doesn't get deleted
    function_.reset( &(*(*system_).function())[index_], dolfin::NoDeleter() );
    oldfunction_.reset( &(*(*system_).oldfunction())[index_], dolfin::NoDeleter() );
    iteratedfunction_.reset( &(*(*system_).iteratedfunction())[index_], dolfin::NoDeleter() );
  }

  // register a pointer to the field as well (first give it a sensible name and label)
  buffer.str(""); buffer << (*system_).name() << "::" << name();
  (*function_).rename(buffer.str(), buffer.str());
  buffer.str(""); buffer << (*system_).name() << "::Old" << name();
  (*oldfunction_).rename(buffer.str(), buffer.str());
  buffer.str(""); buffer << (*system_).name() << "::Iterated" << name();
  (*iteratedfunction_).rename(buffer.str(), buffer.str());

  buffer.str(""); buffer << optionpath() << "/type/rank/boundary_condition";
  int nbcs = Spud::option_count(buffer.str());
  if (nbcs > 0)
  {
    // get the edge id information to set the bcs
    MeshFunction_uint_ptr edgeidmeshfunction = (*(*system_).mesh()).data().mesh_function("EdgeIDs");

    for (uint i = 0; i < nbcs; i++)
    {
      buffer.str(""); buffer << optionpath() << "/type/rank/boundary_condition[" << i << "]";
      bc_fill_(buffer.str(), edgeidmeshfunction);
    }
  }
    
  buffer.str(""); buffer << optionpath() << "/type/rank/initial_condition";
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
      buffer.str(""); buffer << optionpath() << "/type/rank/initial_condition[" << i << "]";
      ic_fill_(buffer.str());
//    }
  }
   
}

void SpudFunctionBucket::bc_fill_(const std::string &optionpath,
                                  const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;
  Spud::OptionError serr;
  
  std::string bcname;
  buffer.str(""); buffer << optionpath << "/name";
  serr = Spud::get_option(buffer.str(), bcname); spud_err(buffer.str(), serr);
  
  std::vector<int> bcids;
  buffer.str(""); buffer << optionpath << "/boundary_ids";
  serr = Spud::get_option(buffer.str(), bcids); spud_err(buffer.str(), serr);
  
  buffer.str(""); buffer << optionpath << "/sub_components";
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsubcomp; i++)
  {
    buffer.str(""); buffer << optionpath << "/sub_components[" << i << "]";
    bc_component_fill_(buffer.str(), bcname, bcids, edgeidmeshfunction);
  }
}

void SpudFunctionBucket::bc_component_fill_(const std::string &optionpath,
                                            const std::string &bcname,
                                            const std::vector<int> &bcids,
                                            const MeshFunction_uint_ptr &edgeidmeshfunction)
{
  std::stringstream buffer;
  std::stringstream namebuffer;
  Spud::OptionError serr;

  buffer.str(""); buffer << optionpath << "/components";
  if (Spud::have_option(buffer.str()))
  {
    // FIXME: tensor support needs to go in here, so another switch between a list and a python function
    //        to return the components!
    std::vector<int> subcompids;
    serr = Spud::get_option(buffer.str(), subcompids); spud_err(buffer.str(), serr);
    
    for (std::vector<int>::const_iterator subcompid = subcompids.begin(); subcompid < subcompids.end(); subcompid++)
    {

       FunctionSpace_ptr subfunctionspace;
       if (contains_subfunctionspace(*subcompid))
       {
         subfunctionspace = fetch_subfunctionspace(*subcompid);
       }
       else
       {
         subfunctionspace.reset( new dolfin::SubSpace(*(functionspace()), *subcompid) );
         register_subfunctionspace(subfunctionspace, *subcompid);
       }
       
       buffer.str(""); buffer << optionpath << "/type::Dirichlet";
       Expression_ptr bcexp = initialize_expression(buffer.str(), size_, shape_);
       
       namebuffer.str(""); namebuffer << bcname << "::" << *subcompid;
       register_bcexpression(bcexp, namebuffer.str());
       
       for (std::vector<int>::const_iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
       {
         BoundaryCondition_ptr bc(new dolfin::DirichletBC(*subfunctionspace, *bcexp, *edgeidmeshfunction, *bcid));
         namebuffer.str(""); namebuffer << bcname << "::" << *subcompid << "::" << *bcid;
         register_bc(bc, namebuffer.str());
       }
       
    }
  }
  else
  {
    buffer.str(""); buffer << optionpath << "/type::Dirichlet";
    Expression_ptr bcexp = initialize_expression(buffer.str(), size_, shape_);
    
    register_bcexpression(bcexp, bcname);
    
    for (std::vector<int>::const_iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
    {
      BoundaryCondition_ptr bc(new dolfin::DirichletBC(*functionspace(), *bcexp, *edgeidmeshfunction, *bcid));
      namebuffer.str(""); namebuffer << bcname << "::" << *bcid;
      register_bc(bc, namebuffer.str());
    }
  }
}

void SpudFunctionBucket::ic_fill_(const std::string &optionpath)
{
  // This will get more complicated when we support looping over regions
  // Loop will probably have to move in here!
  icexpression_ = initialize_expression(optionpath, size_, shape_);
}

void SpudFunctionBucket::initialize_expression_coeff_()
{
  std::stringstream buffer;

  if (type_=="Expression" || type_=="Constant")
  {

    buffer.str(""); buffer << optionpath() << "/type/rank/value";
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
        buffer.str(""); buffer << optionpath() << "/type/rank/value[" << i << "]";

        function_ = initialize_expression(buffer.str(), size_, shape_);
        oldfunction_ = function_;
        iteratedfunction_ = function_;

//      }
    }
  }

}

void SpudFunctionBucket::initialize_function_coeff()
{
  std::stringstream buffer;

  if (type_=="Function")
  {
    if(!(*(*system_).bucket()).contains_coefficientspace(uflsymbol()))
    {
      std::cerr << "Coefficient " << name() << " declared as type::Function ";
      std::cerr << "but its ufl_symbol was not found in any forms or functionals.\n";
      dolfin::error("Unable to allocate coefficient as no matching functionspace found.");
    }

    // only gets handled once we have a functionspace to put it on
    // (either from a solver form or a functional)
    FunctionSpace_ptr coefficientspace = (*(*system_).bucket()).fetch_coefficientspace(uflsymbol());

    function_.reset( new dolfin::Function(*coefficientspace) );

    buffer.str(""); buffer << optionpath() << "/type/rank/value";
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
        buffer.str(""); buffer << optionpath() << "/type/rank/value[" << i << "]";
        Expression_ptr valueexp = initialize_expression(buffer.str(), size_, shape_);

        (*boost::dynamic_pointer_cast< dolfin::Function >(function_)).interpolate(*valueexp);

//      }
    }

    oldfunction_ = function_;
    iteratedfunction_ = function_;

  }

}

void SpudFunctionBucket::functionals_fill_()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  Spud::OptionError serr;
   
  buffer.str(""); buffer << optionpath() << "/output/include_in_diagnostics/functional";
  int nfuncs = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfuncs; i++)
  {
    buffer.str(""); buffer << optionpath() << "/output/include_in_diagnostics/functional[" << i << "]";
    std::string functionaloptionpath = buffer.str();

    std::string functionalname;
    buffer.str(""); buffer << functionaloptionpath << "/name";
    serr = Spud::get_option(buffer.str(), functionalname); spud_err(buffer.str(), serr);

    Form_ptr functional = ufc_fetch_functional((*system_).name(), name(), functionalname, (*system_).mesh());
    register_functional(functional, functionalname, functionaloptionpath);

    // Loop over the functions requested by the form and hook up pointers
    uint ncoeff = (*functional).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*functional).coefficient_name(i);
      // is this a function coefficient?
      if ((*(*system_).bucket()).contains_uflname(uflsymbol))
      {
        // yes, it's a coefficient function... so let's take this opportunity to register
        // its functionspace (if we don't already have it)
        std::string baseuflsymbol = (*(*system_).bucket()).fetch_uflname(uflsymbol);
        if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
        {
          FunctionSpace_ptr coefficientspace;
          coefficientspace = ufc_fetch_coefficientspace((*system_).name(), name(), functionalname, baseuflsymbol, (*system_).mesh());
          (*(*system_).bucket()).register_coefficientspace(coefficientspace, baseuflsymbol);
        }
      }

    }
  }

}

// Register a functional in the function
void SpudFunctionBucket::register_functional(Form_ptr functional, std::string name, std::string optionpath)
{
  // First check if a functional with this name already exists
  Form_it f_it = functionals_.find(name);
  if (f_it != functionals_.end())
  {
    // if it does, issue an error
    dolfin::error("Functional named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    functionals_[name]            = functional;
    functional_optionpaths_[name] = optionpath;
  }
}

// Return a string describing the contents of the spudfunction
const std::string SpudFunctionBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << functionals_str(indent);
  return s.str();
}

// Describe the contents of the functional_optionpaths_ map
const std::string SpudFunctionBucket::functionals_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = functional_optionpaths_.begin(); s_it != functional_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Functional " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

const bool SpudFunctionBucket::include_in_diagnostics() const
{
  return Spud::have_option(optionpath()+"/output/include_in_diagnostics");
}

