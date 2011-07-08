
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

void SpudSystem::fields_fill_(const std::string &optionpath)
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
//    buffer.str(""); buffer << funcbuffer.str() << "/boundary_condition";
//    int nbcs = Spud::option_count(buffer.str());
//    for (uint bci = 0; bci < nbcs; bci++)
//    {
//      std::stringstream bcpathstream;
//      bcpathstream.str(""); bcpathstream << funcbuffer.str() << "/boundary_condition[" << bci << "]";
//      
//      bc_fill_(bcpathstream.str(), funci, dimension, subsysspace, sysname, funcname);
//    }
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
//      GenericFunction_ptr icexpr = init_exp_(icpathstream.str(), dimension);
//      
//      icexprs.insert(std::pair< uint, GenericFunction_ptr >(component, icexpr));
//      
//      component += (*icexpr).value_size();
////     }
//    

}

//void SpudSystem::fill(Bucket_ptr bucket)
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
    fields_fill_(buffer.str()); 
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

