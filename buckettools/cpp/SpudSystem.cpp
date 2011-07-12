
#include "DolfinBoostTypes.h"
#include "SpudSystem.h"
#include "System.h"
#include "SystemsWrapper.h"
#include <dolfin.h>
#include <string>
#include <spud.h>

using namespace buckettools;

SpudSystem::SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh) : optionpath_(optionpath), System(name, mesh)
{
  // Do nothing
}

SpudSystem::~SpudSystem()
{
  // Do nothing
}

void SpudSystem::fields_fill_(const std::string &optionpath, const uint &field_i)
{
  std::stringstream buffer;

  std::string fieldname;
  buffer.str(""); buffer << optionpath << "/name";
  Spud::get_option(buffer.str(), fieldname);

//    buffer.str(""); buffer << sysname << "::" << funcname << "::Function";
//    (*sysfunc)[funci].rename(buffer.str(), buffer.str());
//    
//    buffer.str(""); buffer << sysname << "::" << funcname << "::Residual";
//    (*sysresid)[funci].rename(buffer.str(), buffer.str());
//    
//    FunctionSpace_ptr subsysspace(new dolfin::SubSpace(*sysspace, funci));
//    buffer.str(""); buffer << sysname << "::" << funcname << "::Space";
//    register_functionspace(subsysspace, buffer.str());
//    
  buffer.str(""); buffer << optionpath << "/boundary_condition";
  int nbcs = Spud::option_count(buffer.str());
  for (uint i = 0; i < nbcs; i++)
  {
    buffer << "[" << i << "]";
    bc_fill_(buffer.str(), field_i, dimension, subsysspace, fieldname);
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

void SpudBucket::bc_fill_(const std::string &optionpath, 
                          const int &field_i, 
                          const int dimension, 
                          FunctionSpace_ptr subsysspace,
                          const std::string &fieldname)
{
  std::stringstream buffer;
  
//  MeshFunction_uint_ptr edge_subdomain = fetch_meshfunction("EdgeIDs");
  
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
    std::stringstream subcompbuffer;
    subcompbuffer.str(""); subcompbuffer << optionpath << "/sub_components[" << i << "]";
    
    std::vector<int> subcompids;
    buffer.str(""); buffer << subcompbuffer.str() << "/components";
    if (Spud::have_option(buffer.str()))
    {
      Spud::get_option(buffer.str(), subcompids);
      
      for (std::vector<int>::iterator subcompid = subcompids.begin(); subcompid < subcompids.end(); subcompid++)
      {
        
        FunctionSpace_ptr subsubsysspace( new dolfin::SubSpace(*subsysspace, *subcompid) );
        buffer.str(""); buffer << sysname << "::" << funcname << "::Space_" << *subcompid;
        register_functionspace(subsubsysspace, buffer.str());
        
        buffer.str(""); buffer << subcompbuffer.str() << "/type::Dirichlet";
        GenericFunction_ptr bcexp = initialize_expression(buffer.str(), dimension);
        
        register_bcexp(bcexp);
        
        for (std::vector<int>::iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
        {
          DirichletBC_ptr bc(new dolfin::DirichletBC(*subsubsysspace, *bcexp, *edge_subdomain, *bcid));
          buffer.str(""); buffer << bcname << "_" << *bcid << "_" << *subcompid;
          register_dirichletbc(bc, buffer.str());
        }
        
      }
    }
    else
    {
      buffer.str(""); buffer << subcompbuffer.str() << "/type::Dirichlet";
      GenericFunction_ptr bcexp = initialize_expression(buffer.str(), dimension);
      
      register_bcexp(bcexp);
      
      for (std::vector<int>::iterator bcid = bcids.begin(); bcid < bcids.end(); bcid++)
      {
        DirichletBC_ptr bc(new dolfin::DirichletBC(*subsysspace, *bcexp, *edge_subdomain, *bcid));
        buffer.str(""); buffer << bcname << "_" << *bcid;
        register_dirichletbc(bc, buffer.str());
      }
    }
  }
}

void SpudSystem::fill()
{
  // A buffer to put option paths in
  std::stringstream buffer;

  functionspace_ = buckettools::fetch_functionspace(name(), mesh_);

  function_.reset( new dolfin::Function(*functionspace_) );
  oldfunction_.reset( new dolfin::Function(*functionspace_) );
  iteratedfunction_.reset( new dolfin::Function(*functionspace_) );

//  std::map< uint, GenericFunction_ptr > icexprs;
//  uint component = 0;
//  
  buffer.str("");  buffer << optionpath() << "/field";
  int nfields = Spud::option_count(buffer.str());
  for (uint i = 0; i < nfields; i++)
  {
    buffer << "[" << i << "]";
    fields_fill_(buffer.str(), i); 
  }
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

