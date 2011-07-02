
#include "SystemBucket.h"
#include "Bucket.h"
#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

SystemBucket::SystemBucket()
{
  // Do nothing
}

SystemBucket::~SystemBucket()
{
  // Do nothing
}

void SystemBucket::register_mesh(Mesh_ptr mesh, std::string name)
{
  meshes_.insert(std::pair<std::string, Mesh_ptr>(name, mesh));
}

void SystemBucket::register_meshfunction(MeshFunction_uint_ptr meshfunction, std::string name)
{
  meshfunctions_.insert(std::pair<std::string, MeshFunction_uint_ptr>(name, meshfunction));
}

void SystemBucket::register_functionspace(FunctionSpace_ptr functionspace, std::string name)
{
  functionspaces_.insert(std::pair<std::string, FunctionSpace_ptr>(name, functionspace));
}

void SystemBucket::register_dirichletbc(DirichletBC_ptr dirichletbc, std::string name)
{
  dirichletbcs_.insert(std::pair<std::string, DirichletBC_ptr>(name, dirichletbc));
}

void SystemBucket::register_detector(Detectors_ptr detector, std::string name)
{
  detectors_.insert(std::pair<std::string, Detectors_ptr>(name, detector));
}

void SystemBucket::register_function(Function_ptr function, std::string name)
{
  functions_.insert(std::pair<std::string, Function_ptr>(name, function));
}

void SystemBucket::register_bcexp(GenericFunction_ptr bcexp)
{
  bcexps_.push_back(bcexp);
}

Mesh_ptr SystemBucket::fetch_mesh(const std::string name)
{
  std::map< std::string, Mesh_ptr >::iterator it;
  it = meshes_.find(name);
  return (*it).second;
}

MeshFunction_uint_ptr SystemBucket::fetch_meshfunction(const std::string name)
{
  std::map< std::string, MeshFunction_uint_ptr >::iterator it;
  it = meshfunctions_.find(name);
  return (*it).second;
}

FunctionSpace_ptr SystemBucket::fetch_functionspace(const std::string name)
{
  std::map< std::string, FunctionSpace_ptr >::iterator it;
  it = functionspaces_.find(name);
  return (*it).second;
}

Detectors_ptr SystemBucket::fetch_detector(const std::string name)
{
  std::map< std::string, Detectors_ptr >::iterator it;
  it = detectors_.find(name);
  return (*it).second;
}

Function_ptr SystemBucket::fetch_function(const std::string name)
{
  std::map< std::string, Function_ptr >::iterator it;
  it = functions_.find(name);
  return (*it).second;
}

std::map< std::string, Detectors_ptr >::iterator SystemBucket::detectors_begin()
{
  return detectors_.begin();
}

std::map< std::string, Detectors_ptr >::const_iterator SystemBucket::detectors_begin() const
{
  return detectors_.begin();
}

std::map< std::string, Detectors_ptr >::iterator SystemBucket::detectors_end()
{
  return detectors_.end();
}

std::map< std::string, Detectors_ptr >::const_iterator SystemBucket::detectors_end() const
{
  return detectors_.end();
}

std::map< std::string, DirichletBC_ptr >::iterator SystemBucket::dirichletbcs_begin()
{
  return dirichletbcs_.begin();
}

std::map< std::string, DirichletBC_ptr >::const_iterator SystemBucket::dirichletbcs_begin() const
{
  return dirichletbcs_.begin();
}

std::map< std::string, DirichletBC_ptr >::iterator SystemBucket::dirichletbcs_end()
{
  return dirichletbcs_.end();
}

std::map< std::string, DirichletBC_ptr >::const_iterator SystemBucket::dirichletbcs_end() const
{
  return dirichletbcs_.end();
}

std::map< std::string, Function_ptr >::iterator SystemBucket::functions_begin()
{
  return functions_.begin();
}

std::map< std::string, Function_ptr >::const_iterator SystemBucket::functions_begin() const
{
  return functions_.begin();
}

std::map< std::string, Function_ptr >::iterator SystemBucket::functions_end()
{
  return functions_.end();
}

std::map< std::string, Function_ptr >::const_iterator SystemBucket::functions_end() const
{
  return functions_.end();
}

void SystemBucket::clean_()
{
  
  meshes_.clear();
  meshfunctions_.clear();
  detectors_.clear();
  functions_.clear();
  
}
