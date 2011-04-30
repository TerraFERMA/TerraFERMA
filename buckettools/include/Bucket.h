
#ifndef __BUCKET_H
#define __BUCKET_H

#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>

namespace buckettools
{
  
  typedef boost::shared_ptr< dolfin::GenericFunction > GenericFunction_ptr;
  typedef boost::shared_ptr< dolfin::Mesh > Mesh_ptr;
  typedef boost::shared_ptr< dolfin::MeshFunction< dolfin::uint > > MeshFunction_uint_ptr;
  typedef boost::shared_ptr< dolfin::FunctionSpace > FunctionSpace_ptr;
  typedef boost::shared_ptr< dolfin::DirichletBC > DirichletBC_ptr;
  
  class Bucket
  {
  private:
    
    std::map<std::string, Mesh_ptr> meshes_;
    
    std::map< std::string, MeshFunction_uint_ptr > meshfunctions_;
    
    std::map< std::string, FunctionSpace_ptr > functionspaces_;
    
    std::map< std::string, DirichletBC_ptr > dirichletbcs_;
    
    std::map< std::string, GenericFunction_ptr > functions_;
    
    std::map< std::string, Detectors_ptr > detectors_;
    
    
    void clean_();
    
  public:
    
    Bucket();
    
    ~Bucket();
    
    void register_mesh(Mesh_ptr mesh, std::string name);
    
    void register_meshfunction(MeshFunction_uint_ptr meshfunction, std::string name);
    
    void register_functionspace(FunctionSpace_ptr functionspace, std::string name);
    
    void register_dirichletbc(DirichletBC_ptr dirichletbc, std::string name);
    
    void register_detector(Detectors_ptr detector, std::string name);
    
    void register_function(GenericFunction_ptr function, std::string name);
    
    Mesh_ptr fetch_mesh(const std::string name);
    
    MeshFunction_uint_ptr fetch_meshfunction(const std::string name);
    
    FunctionSpace_ptr fetch_functionspace(const std::string name);
    
    GenericFunction_ptr fetch_function(const std::string name);
    
    Detectors_ptr fetch_detector(const std::string name);
    
    std::map< std::string, Detectors_ptr >::iterator detectors_begin();
    
    std::map< std::string, Detectors_ptr >::const_iterator detectors_begin() const;
    
    std::map< std::string, Detectors_ptr >::iterator detectors_end();
    
    std::map< std::string, Detectors_ptr >::const_iterator detectors_end() const;
    
    std::map< std::string, GenericFunction_ptr >::iterator functions_begin();
    
    std::map< std::string, GenericFunction_ptr >::const_iterator functions_begin() const;
    
    std::map< std::string, GenericFunction_ptr >::iterator functions_end();
    
    std::map< std::string, GenericFunction_ptr >::const_iterator functions_end() const;
    
  };
}
#endif