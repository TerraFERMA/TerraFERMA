
#include "DolfinBoostTypes.h"
#include "SpudSystem.h"
#include "System.h"
#include "SystemsWrapper.h"
#include "SpudBase.h"
#include <dolfin.h>
#include <string>
#include <spud.h>

using namespace buckettools;

// Specific constructor
SpudSystem::SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh) : optionpath_(optionpath), System(name, mesh)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudSystem::~SpudSystem()
{
  // Do nothing
}

// Fill the system using spud and assuming a buckettools schema structure
void SpudSystem::fill(const uint &dimension)
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;

  // Here's where the automatically generated magic happens... this is fetching
  // the functionspace from ufl
  functionspace_ = buckettools::fetch_functionspace(name(), mesh_);

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

  buffer.str("");  buffer << optionpath() << "/field";
  int nfields = Spud::option_count(buffer.str());
  // Loop over the fields (which are subfunctions of this functionspace)
  // and register them in the system
  for (uint i = 0; i < nfields; i++)
  {
    buffer << "[" << i << "]";
    fields_fill_(buffer.str(), i, nfields, dimension); 
  }

//  
//  std::map< uint, GenericFunction_ptr > icexprs;
//  uint component = 0;
//  
//  InitialConditionExpression sysicexpr(component, icexprs);
//  (*sysfunc).interpolate(sysicexpr);
//  for(std::map< std::string, DirichletBC_ptr >::iterator
//            bc = dirichletbcs_begin(); 
//            bc != dirichletbcs_end(); bc++)
//  {
//    (*((*bc).second)).apply((*sysfunc).vector());
//  }
//  
}

// Fill out the information regarding the each subfunction (or field)
void SpudSystem::fields_fill_(const std::string &optionpath, 
                              const uint &field_i, 
                              const uint &nfields, 
                              const uint &dimension)
{
  // A buffer to put option paths in
  std::stringstream buffer;

  // What is the field size (if it's a vector)
  // Would it be possible to get this from the subfunctionspace below?
  int size;
  buffer.str(""); buffer << optionpath << "/type/rank/element/size";
  Spud::get_option(buffer.str(), size, dimension);

  // What is the field shape (if it's a tensor)
  // Would it be possible to get this from the subfunctionspace below?
  std::vector< int > shape;
  std::vector< int > default_shape(2, dimension);
  buffer.str(""); buffer << optionpath << "/type/rank/element/shape";
  Spud::get_option(buffer.str(), shape, default_shape);

  // What is this field called?
  std::string fieldname;
  buffer.str(""); buffer << optionpath << "/name";
  Spud::get_option(buffer.str(), fieldname);

  // Is this a mixed functionspace or not?
  FunctionSpace_ptr subfunctionspace;
  if (nfields == 1)
  {
    // no... the subfunctionspace for this field is identical to the system's
    // luckily for us these are just pointers so grab a reference to it
    subfunctionspace = functionspace_;
  }
  else
  {
    // yes... use DOLFIN to extract that subspace so we can declare things on it (ics, bcs etc.)
    subfunctionspace.reset( new dolfin::SubSpace(*functionspace_, field_i) );
  }
  // register a pointer to the subfunctionspace in the system so it remains in memory for other
  // objects that take a reference
  register_subfunctionspace(subfunctionspace, fieldname);

  buffer.str(""); buffer << optionpath << "/type/rank/boundary_condition";
  int nbcs = Spud::option_count(buffer.str());
  if (nbcs > 0)
  {
    // get the edge id information to set the bcs
    MeshFunction_uint_ptr edgeidmeshfunction = (*mesh_).data().mesh_function("EdgeIDs");

    for (uint i = 0; i < nbcs; i++)
    {
      buffer << "[" << i << "]";
      bc_fill_(buffer.str(), fieldname, size, shape, subfunctionspace, edgeidmeshfunction);
    }
  }
//    
//    buffer.str(""); buffer << funcbuffer.str() << "/initial_condition";
//    int nics = Spud::option_count(buffer.str());
//    if (nics > 1)
//    {
//      dolfin::error("Haven't thought about ics over regions.");
//    }
////     for (uint ici = 0; ici < nics; ici++)
////     {
//      uint ici = 0;
//      
//      std::stringstream icpathstream;
//      icpathstream.str(""); icpathstream << funcbuffer.str() << "/initial_condition[" << ici << "]";
//      
//      GenericFunction_ptr icexpr = initialize_expression(icpathstream.str(), dimension);
//      
//      icexprs.insert(std::pair< uint, GenericFunction_ptr >(component, icexpr));
//      
//      component += (*icexpr).value_size();
////     }
//    

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
    buffer << "[" << i << "]";
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
    
    for (std::vector<int>::iterator subcompid = subcompids.begin(); subcompid < subcompids.end(); subcompid++)
    {
//        
//        FunctionSpace_ptr subsubsysspace( new dolfin::SubSpace(*subsysspace, *subcompid) );
//        buffer.str(""); buffer << sysname << "::" << funcname << "::Space_" << *subcompid;
//        register_functionspace(subsubsysspace, buffer.str());
//        
//        buffer.str(""); buffer << subcompbuffer.str() << "/type::Dirichlet";
//        GenericFunction_ptr bcexp = initialize_expression(buffer.str(), dimension);
//        
//        register_bcexp(bcexp);
//        
//        for (std::vector<int>::iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
//        {
//          DirichletBC_ptr bc(new dolfin::DirichletBC(*subsubsysspace, *bcexp, *edge_subdomain, *bcid));
//          buffer.str(""); buffer << bcname << "_" << *bcid << "_" << *subcompid;
//          register_dirichletbc(bc, buffer.str());
//        }
//        
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

