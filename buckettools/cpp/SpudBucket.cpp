
#include "SpudBucket.h"
#include "Detectors.h"
#include "PointDetectors.h"
#include "PythonDetectors.h"
#include "PythonExpression.h"
#include "InitialConditionExpression.h"
#include "System.h"
#include <dolfin.h>
#include <string>
#include <spud.h>

using namespace buckettools;

SpudBucket::SpudBucket() : Bucket()
{
  // Do nothing
}

SpudBucket::~SpudBucket()
{
  // Do nothing
}

void SpudBucket::fill()
{
  
  meshes_fill_();
  
  int nsystems = Spud::option_count("/system");
  for (uint i = 0; i<nsystems; i++)
  {
    system_fill_(i);
  }
  
  detectors_fill_();
  
}

void SpudBucket::meshes_fill_()
{
  
  std::string filename;
  
  Spud::get_option("/geometry/mesh/file", filename);
  
  Mesh_ptr mesh(new dolfin::Mesh(filename));
  (*mesh).init();
  register_mesh(mesh, "Mesh");
  
  Spud::get_option("/geometry/mesh/edges_file", filename);
  
  MeshFunction_uint_ptr meshfunc(new dolfin::MeshFunction<dolfin::uint>(*mesh, filename));
  register_meshfunction(meshfunc, "EdgeIDs");
  
  Spud::get_option("/geometry/mesh/cells_file", filename);
  
  meshfunc.reset(new dolfin::MeshFunction<dolfin::uint>(*mesh, filename));
  register_meshfunction(meshfunc, "CellIDs");
  
}

void SpudBucket::system_fill_(const uint &sysindex)
{
  std::stringstream buffer;
  
  std::string sysname;
  buffer.str(""); buffer << "/system[" << sysindex << "]/name";
  Spud::get_option(buffer.str(), sysname);
  
  int dimension;
  Spud::get_option("/geometry/dimension", dimension);
  
  Mesh_ptr mesh = fetch_mesh("Mesh");
  
  FunctionSpace_ptr sysspace(new System::FunctionSpace(mesh));
  buffer.str(""); buffer << sysname << "::Space";
  register_functionspace(sysspace, buffer.str());
  
  Function_ptr sysfunc(new dolfin::Function(*sysspace));
  buffer.str(""); buffer << sysname << "::Function";
  (*sysfunc).rename(buffer.str(), buffer.str());
  register_function(sysfunc, buffer.str());
  Function_ptr sysresid(new dolfin::Function(*sysspace));
  buffer.str(""); buffer << sysname << "::Residual";
  (*sysresid).rename(buffer.str(), buffer.str());
  register_function(sysresid, buffer.str());
  
  std::map< uint, GenericFunction_ptr > icexprs;
  uint component = 0;
  
  buffer.str("");  buffer << "/system[" << sysindex << "]/function";
  int nfunctions = Spud::option_count(buffer.str());
  for (uint funci = 0; funci < nfunctions; funci++)
  {
    std::stringstream funcbuffer;
    funcbuffer.str(""); funcbuffer << "/system[" << sysindex << "]/function[" << funci << "]";
    
    std::string funcname;
    buffer.str(""); buffer << funcbuffer.str() << "/name";
    Spud::get_option(buffer.str(), funcname);
    
    buffer.str(""); buffer << sysname << "::" << funcname << "::Function";
    (*sysfunc)[funci].rename(buffer.str(), buffer.str());
    
    buffer.str(""); buffer << sysname << "::" << funcname << "::Residual";
    (*sysresid)[funci].rename(buffer.str(), buffer.str());
    
    FunctionSpace_ptr subsysspace(new dolfin::SubSpace(*sysspace, funci));
    buffer.str(""); buffer << sysname << "::" << funcname << "::Space";
    register_functionspace(subsysspace, buffer.str());
    
    buffer.str(""); buffer << funcbuffer.str() << "/boundary_condition";
    int nbcs = Spud::option_count(buffer.str());
    for (uint bci = 0; bci < nbcs; bci++)
    {
      std::stringstream bcpathstream;
      bcpathstream.str(""); bcpathstream << funcbuffer.str() << "/boundary_condition[" << bci << "]";
      
      bc_fill_(bcpathstream.str(), funci, dimension, subsysspace, sysname, funcname);
    }
    
    buffer.str(""); buffer << funcbuffer.str() << "/initial_condition";
    int nics = Spud::option_count(buffer.str());
    if (nics > 1)
    {
      dolfin::error("Haven't thought about ics over regions.");
    }
//     for (uint ici = 0; ici < nics; ici++)
//     {
      uint ici = 0;
      
      std::stringstream icpathstream;
      icpathstream.str(""); icpathstream << funcbuffer.str() << "/initial_condition[" << ici << "]";
      
      GenericFunction_ptr icexpr = init_exp_(icpathstream.str(), dimension);
      
      icexprs.insert(std::pair< uint, GenericFunction_ptr >(component, icexpr));
      
      component += (*icexpr).value_size();
//     }
    
  }
  
  InitialConditionExpression sysicexpr(component, icexprs);
  (*sysfunc).interpolate(sysicexpr);
  for(std::map< std::string, DirichletBC_ptr >::iterator
            bc = dirichletbcs_begin(); 
            bc != dirichletbcs_end(); bc++)
  {
    (*((*bc).second)).apply((*sysfunc).vector());
  }
  
}

void SpudBucket::bc_fill_(const std::string bcpath, 
                          const int funci, 
                          const int dimension, 
                          FunctionSpace_ptr subsysspace,
                          const std::string sysname,
                          const std::string funcname)
{
  std::stringstream buffer;
  
  MeshFunction_uint_ptr edge_subdomain = fetch_meshfunction("EdgeIDs");
  
  std::string bcname;
  buffer.str(""); buffer << bcpath << "/name";
  Spud::get_option(buffer.str(), bcname);
  
  std::vector<int> bcids;
  buffer.str(""); buffer << bcpath << "/boundary_ids";
  Spud::get_option(buffer.str(), bcids);
  
  buffer.str(""); buffer << bcpath << "/sub_components";
  int nsubcomp = Spud::option_count(buffer.str());
  for (uint subcompi = 0; subcompi < nsubcomp; subcompi++)
  {
    std::stringstream subcompbuffer;
    subcompbuffer.str(""); subcompbuffer << bcpath << "/sub_components[" << subcompi << "]";
    
    std::vector<int> subcompids;
    buffer.str(""); buffer << subcompbuffer.str() << "/components";
    if(Spud::have_option(buffer.str()))
    {
      Spud::get_option(buffer.str(), subcompids);
      
      for (std::vector<int>::iterator subcompid = subcompids.begin(); subcompid < subcompids.end(); subcompid++)
      {
        
        FunctionSpace_ptr subsubsysspace(new dolfin::SubSpace(*subsysspace, *subcompid));
        buffer.str(""); buffer << sysname << "::" << funcname << "::Space_" << *subcompid;
        register_functionspace(subsubsysspace, buffer.str());
        
        buffer.str(""); buffer << subcompbuffer.str() << "/type::Dirichlet";
        GenericFunction_ptr bcexp = init_exp_(buffer.str(), dimension);
        
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
      GenericFunction_ptr bcexp = init_exp_(buffer.str(), dimension);
      
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

GenericFunction_ptr SpudBucket::init_exp_(const std::string path, const int dimension)
{
  GenericFunction_ptr bcexp(new dolfin::Expression());
  
  std::stringstream constbuffer, pybuffer;
  constbuffer.str(""); constbuffer << path << "/constant";
  pybuffer.str(""); pybuffer << path << "/python";
  
  if(Spud::have_option(constbuffer.str()))
  {
    int rank;
    Spud::get_option_rank(constbuffer.str(), rank);
    if(rank==0)
    {
      double value;
      Spud::get_option(constbuffer.str(), value);
      bcexp.reset(new dolfin::Constant(value));
    }
    else if (rank==1)
    {
      std::vector<double> values;
      Spud::get_option(constbuffer.str(), values);
      bcexp.reset(new dolfin::Constant(values));
    }
    else
    {
      dolfin::error("Don't know how to deal with rank > 1 yet.");
    }
  } 
  else if(Spud::have_option(pybuffer.str()))
  {
    std::string pyfunction;
    Spud::get_option(pybuffer.str(), pyfunction);
    
    // rank of a python function isn't in the default spud base language
    // so have added it... but it comes out as a string of course!
    std::stringstream buffer;
    std::string rankstring; // bit of a hack
    buffer.str(""); buffer << pybuffer.str() << "/rank";
    Spud::get_option(buffer.str(), rankstring);
    
    int rank;
    rank = atoi(rankstring.c_str());
    if(rank==0)
    {
      bcexp.reset(new buckettools::PythonExpression(pyfunction));
    }
    else if (rank==1)
    {
      bcexp.reset(new buckettools::PythonExpression(dimension, pyfunction));
    }
    else
    {
      dolfin::error("Don't know how to deal with rank > 1 yet.");
    }
  }
  else
  {
    dolfin::error("Unknown way of specifying bc expression.");
  }
  
  return bcexp;
  
}

void SpudBucket::detectors_fill_()
{
  int dimension;
  Spud::get_option("/geometry/dimension", dimension);

  Detectors_ptr det(new Detectors());
  std::stringstream buffer;
  
  int ndets;
  ndets = Spud::option_count("/io/detectors/point");
  for (uint i=0; i<ndets; i++)
  {
    std::vector<double> point;
    std::string detname;
    
    buffer.str(""); buffer << "/io/detectors/point[" << i << "]/name";
    Spud::get_option(buffer.str(), detname);
    buffer.str(""); buffer << "/io/detectors/point[" << i << "]";
    Spud::get_option(buffer.str(), point);
    
    det.reset(new PointDetectors(&point, detname));
    register_detector(det, detname);
  }
  
  ndets = Spud::option_count("/io/detectors/array");
  for (uint i=0; i<ndets; i++)
  {
    int no_det;
    std::string detname, function;
    
    buffer.str(""); buffer << "/io/detectors/array[" << i << "]/number_of_detectors";
    Spud::get_option(buffer.str(), no_det);
    buffer.str(""); buffer << "/io/detectors/array[" << i << "]/name";
    Spud::get_option(buffer.str(), detname);
    buffer.str(""); buffer << "/io/detectors/array[" << i << "]/python";
    Spud::get_option(buffer.str(), function);
    
    det.reset(new PythonDetectors(no_det, dimension, function, detname));
    register_detector(det, detname);
  }  
  
}
