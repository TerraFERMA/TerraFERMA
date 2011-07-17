
#include "Bucket.h"
// #include "GenericDetectors.h"
// #include "PythonDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
Bucket::Bucket()
{
  // Do nothing
}

// Specific constructor
Bucket::Bucket(std::string name) : name_(name)
{
  // Do nothing
}

// Default destructor
Bucket::~Bucket()
{
  empty_();
}

// Return a string describing the contents of the bucket
std::string Bucket::str() const
{
  std::stringstream s;
  int indent = 1;
  s << "Bucket " << name() << std::endl;
  s << meshes_str(indent);
  s << systems_str(indent);
  return s.str();
}

// Return a string describing the contents of meshes_
std::string Bucket::meshes_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Mesh_const_it m_it = meshes_begin(); m_it != meshes_end(); m_it++ )
  {
    s << indentation << "Mesh " << (*m_it).first  << std::endl;
  }
  return s.str();
}

// Return a string describing the contents of systems_
std::string Bucket::systems_str(int indent) const
{
  std::stringstream s;
  for ( System_const_it s_it = systems_begin(); s_it != systems_end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
  }
  return s.str();
}

// Register a pointer to a DOLFIN mesh
void Bucket::register_mesh(Mesh_ptr mesh, std::string name)
{
  // First check if a mesh with this name already exists
  Mesh_it m_it = meshes_.find(name);
  if (m_it != meshes_end())
  {
    // if it does, issue an error
    dolfin::error("Mesh named \"%s\" already exists in bucket.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    meshes_[name] = mesh;
  }
}

// Fetch a pointer to a DOLFIN mesh
Mesh_ptr Bucket::fetch_mesh(const std::string name)
{
  // Check if the mesh exists in the bucket
  Mesh_it m_it = meshes_.find(name);
  if (m_it == meshes_end())
  {
    // if it doesn't then throw an error
    dolfin::error("Mesh named \"%s\" does not exist in bucket.", name.c_str());
  }
  else
  {
    // if it does then return the pointer to the mesh
    return (*m_it).second;
  }
}

// Public iterator access
Mesh_it Bucket::meshes_begin()
{
  return meshes_.begin();
}

// Public iterator access
Mesh_const_it Bucket::meshes_begin() const
{
  return meshes_.begin();
}

// Public iterator access
Mesh_it Bucket::meshes_end()
{
  return meshes_.end();
}

// Public iterator access
Mesh_const_it Bucket::meshes_end() const
{
  return meshes_.end();
}

// Register a pointer to a DOLFIN mesh
void Bucket::register_system(System_ptr system, std::string name)
{
  // First check if a system with this name already exists
  System_it s_it = systems_.find(name);
  if (s_it != systems_end())
  {
    // if it does, issue an error
    dolfin::error("System named \"%s\" already exists in bucket.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    systems_[name] = system;
  }
}

// Fetch a pointer to a DOLFIN mesh
System_ptr Bucket::fetch_system(const std::string name)
{
  // Check if the system exists in the bucket
  System_it s_it = systems_.find(name);
  if (s_it == systems_end())
  {
    // if it doesn't then throw an error
    dolfin::error("System named \"%s\" does not exist in bucket.", name.c_str());
  }
  else
  {
    // if it does then return the pointer to the mesh
    return (*s_it).second;
  }
}

// Public iterator access
System_it Bucket::systems_begin()
{
  return systems_.begin();
}

// Public iterator access
System_const_it Bucket::systems_begin() const
{
  return systems_.begin();
}

// Public iterator access
System_it Bucket::systems_end()
{
  return systems_.end();
}

// Public iterator access
System_const_it Bucket::systems_end() const
{
  return systems_.end();
}

// Empty the bucket
void Bucket::empty_()
{
  meshes_.clear();
  systems_.clear();
//  detectors_.clear();
}

//// Create a system bucket in this bucket
//void Bucket::register_system(std::string name, std::string option_path)
//{
//  // First check if a system with this name already exists
//  SystemBucket_it s_it = systembuckets_.find(name);
//  if(s_it != systembuckets_.end())
//  {
//    // if it does, issue an error
//    dolfin::error("System named \"%s\" already exists in bucket.", name.c_str());
//  }
//  else
//  {
//    // if not then create it and insert it into the maps
//    SystemBucket_ptr system(new SystemBucket(name, option_path));
//    systembuckets_[name] = system;
//    systembucket_optionpaths_[name] = option_path;
//  }
//}
//
//// Register a detector set in the bucket
//void Bucket::register_detector(GenericDetectors_ptr detector, std::string name, std::string option_path)
//{
//  // First check if a detector set with this name already exists
//  GenericDetectors_it d_it = detectors_.find(name);
//  if (d_it != detectors_.end())
//  {
//    // if it does, issue an error
//    dolfin::error("Detector set named \"%s\" already exists in bucket.", name.c_str());
//  }
//  else
//  {
//    // if not then insert it into the maps
//    detectors_[name] = detector;
//    detector_optionpaths_[name] = option_path;
//  }
//}
//
//void Bucket::register_functionspace(FunctionSpace_ptr functionspace, std::string name)
//{
//  functionspaces_.insert(std::pair<std::string, FunctionSpace_ptr>(name, functionspace));
//}
//
//void Bucket::register_dirichletbc(DirichletBC_ptr dirichletbc, std::string name)
//{
//  dirichletbcs_.insert(std::pair<std::string, DirichletBC_ptr>(name, dirichletbc));
//}
//
//void Bucket::register_function(Function_ptr function, std::string name)
//{
//  functions_.insert(std::pair<std::string, Function_ptr>(name, function));
//}
//
//void Bucket::register_bcexp(GenericFunction_ptr bcexp)
//{
//  bcexps_.push_back(bcexp);
//}

//FunctionSpace_ptr Bucket::fetch_functionspace(const std::string name)
//{
//  std::map< std::string, FunctionSpace_ptr >::iterator it;
//  it = functionspaces_.find(name);
//  return (*it).second;
//}
//
//GenericDetectors_ptr Bucket::fetch_detector(const std::string name)
//{
//  std::map< std::string, GenericDetectors_ptr >::iterator it;
//  it = detectors_.find(name);
//  return (*it).second;
//}
//
//Function_ptr Bucket::fetch_function(const std::string name)
//{
//  std::map< std::string, Function_ptr >::iterator it;
//  it = functions_.find(name);
//  return (*it).second;
//}
//
//
//std::map< std::string, GenericDetectors_ptr >::iterator Bucket::detectors_begin()
//{
//  return detectors_.begin();
//}
//
//std::map< std::string, GenericDetectors_ptr >::const_iterator Bucket::detectors_begin() const
//{
//  return detectors_.begin();
//}
//
//std::map< std::string, GenericDetectors_ptr >::iterator Bucket::detectors_end()
//{
//  return detectors_.end();
//}
//
//std::map< std::string, GenericDetectors_ptr >::const_iterator Bucket::detectors_end() const
//{
//  return detectors_.end();
//}
//
//std::map< std::string, DirichletBC_ptr >::iterator Bucket::dirichletbcs_begin()
//{
//  return dirichletbcs_.begin();
//}
//
//std::map< std::string, DirichletBC_ptr >::const_iterator Bucket::dirichletbcs_begin() const
//{
//  return dirichletbcs_.begin();
//}
//
//std::map< std::string, DirichletBC_ptr >::iterator Bucket::dirichletbcs_end()
//{
//  return dirichletbcs_.end();
//}
//
//std::map< std::string, DirichletBC_ptr >::const_iterator Bucket::dirichletbcs_end() const
//{
//  return dirichletbcs_.end();
//}
//
//std::map< std::string, Function_ptr >::iterator Bucket::functions_begin()
//{
//  return functions_.begin();
//}
//
//std::map< std::string, Function_ptr >::const_iterator Bucket::functions_begin() const
//{
//  return functions_.begin();
//}
//
//std::map< std::string, Function_ptr >::iterator Bucket::functions_end()
//{
//  return functions_.end();
//}
//
//std::map< std::string, Function_ptr >::const_iterator Bucket::functions_end() const
//{
//  return functions_.end();
//}

