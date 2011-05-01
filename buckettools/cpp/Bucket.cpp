
#include "Bucket.h"
#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

Bucket::Bucket()
{
  // Do nothing
}

Bucket::~Bucket()
{
  // Do nothing
}

void Bucket::register_mesh(Mesh_ptr mesh, std::string name)
{
  meshes_.insert(std::pair<std::string, Mesh_ptr>(name, mesh));
}

void Bucket::register_meshfunction(MeshFunction_uint_ptr meshfunction, std::string name)
{
  meshfunctions_.insert(std::pair<std::string, MeshFunction_uint_ptr>(name, meshfunction));
}

void Bucket::register_functionspace(FunctionSpace_ptr functionspace, std::string name)
{
  functionspaces_.insert(std::pair<std::string, FunctionSpace_ptr>(name, functionspace));
}

void Bucket::register_dirichletbc(DirichletBC_ptr dirichletbc, std::string name)
{
  dirichletbcs_.insert(std::pair<std::string, DirichletBC_ptr>(name, dirichletbc));
}

void Bucket::register_detector(Detectors_ptr detector, std::string name)
{
  detectors_.insert(std::pair<std::string, Detectors_ptr>(name, detector));
}

void Bucket::register_function(Function_ptr function, std::string name)
{
  functions_.insert(std::pair<std::string, Function_ptr>(name, function));
}

void Bucket::register_bcexp(GenericFunction_ptr bcexp)
{
  bcexps_.push_back(bcexp);
}

Mesh_ptr Bucket::fetch_mesh(const std::string name)
{
  std::map< std::string, Mesh_ptr >::iterator it;
  it = meshes_.find(name);
  return (*it).second;
}

MeshFunction_uint_ptr Bucket::fetch_meshfunction(const std::string name)
{
  std::map< std::string, MeshFunction_uint_ptr >::iterator it;
  it = meshfunctions_.find(name);
  return (*it).second;
}

FunctionSpace_ptr Bucket::fetch_functionspace(const std::string name)
{
  std::map< std::string, FunctionSpace_ptr >::iterator it;
  it = functionspaces_.find(name);
  return (*it).second;
}

Detectors_ptr Bucket::fetch_detector(const std::string name)
{
  std::map< std::string, Detectors_ptr >::iterator it;
  it = detectors_.find(name);
  return (*it).second;
}

Function_ptr Bucket::fetch_function(const std::string name)
{
  std::map< std::string, Function_ptr >::iterator it;
  it = functions_.find(name);
  return (*it).second;
}

std::map< std::string, Detectors_ptr >::iterator Bucket::detectors_begin()
{
  return detectors_.begin();
}

std::map< std::string, Detectors_ptr >::const_iterator Bucket::detectors_begin() const
{
  return detectors_.begin();
}

std::map< std::string, Detectors_ptr >::iterator Bucket::detectors_end()
{
  return detectors_.end();
}

std::map< std::string, Detectors_ptr >::const_iterator Bucket::detectors_end() const
{
  return detectors_.end();
}

std::map< std::string, DirichletBC_ptr >::iterator Bucket::dirichletbcs_begin()
{
  return dirichletbcs_.begin();
}

std::map< std::string, DirichletBC_ptr >::const_iterator Bucket::dirichletbcs_begin() const
{
  return dirichletbcs_.begin();
}

std::map< std::string, DirichletBC_ptr >::iterator Bucket::dirichletbcs_end()
{
  return dirichletbcs_.end();
}

std::map< std::string, DirichletBC_ptr >::const_iterator Bucket::dirichletbcs_end() const
{
  return dirichletbcs_.end();
}

std::map< std::string, Function_ptr >::iterator Bucket::functions_begin()
{
  return functions_.begin();
}

std::map< std::string, Function_ptr >::const_iterator Bucket::functions_begin() const
{
  return functions_.begin();
}

std::map< std::string, Function_ptr >::iterator Bucket::functions_end()
{
  return functions_.end();
}

std::map< std::string, Function_ptr >::const_iterator Bucket::functions_end() const
{
  return functions_.end();
}

void Bucket::clean_()
{
  
  meshes_.clear();
  meshfunctions_.clear();
  detectors_.clear();
  functions_.clear();
  
}
