
#include "Bucket.h"
// #include "GenericDetectors.h"
// #include "PythonDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
Bucket::Bucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
Bucket::Bucket(std::string name) : name_(name)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
Bucket::~Bucket()
{
  empty_();                                                          // empty the data structures (unnecessary?)
}

//*******************************************************************|************************************************************//
// run the model
//*******************************************************************|************************************************************//
void Bucket::run()
{
  dolfin::error("Failed to find virtual function run.");             // require a function to be defined in a derived class
}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling solve on each of them
//*******************************************************************|************************************************************//
void Bucket::solve()
{
  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                              s_it != orderedsystems_end(); s_it++)
  {
    (*(*s_it).second).solve();
  }
}

//*******************************************************************|************************************************************//
// loop over the ordered systems in the bucket, calling update on each of them
//*******************************************************************|************************************************************//
void Bucket::update()
{
  for (int_SystemBucket_const_it s_it = orderedsystems_begin(); 
                             s_it != orderedsystems_end(); s_it++)
  {
    (*(*s_it).second).update();
  }
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a dolfin mesh in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_mesh(Mesh_ptr mesh, const std::string &name)
{
  Mesh_it m_it = meshes_.find(name);                                 // check if a mesh with this name already exists
  if (m_it != meshes_end())
  {
    dolfin::error("Mesh named \"%s\" already exists in bucket.",     // if it does, issue an error
                                                    name.c_str());
  }
  else
  {
    meshes_[name] = mesh;                                            // if not, insert it into the meshes_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a dolfin mesh in the bucket data maps
//*******************************************************************|************************************************************//
Mesh_ptr Bucket::fetch_mesh(const std::string &name)
{
  Mesh_it m_it = meshes_.find(name);                                 // check if this mesh exists in the meshes_ map
  if (m_it == meshes_end())
  {
    dolfin::error("Mesh named \"%s\" does not exist in bucket.",     // if it doesn't, issue an error
                                                    name.c_str());
  }
  else
  {
    return (*m_it).second;                                           // if it does, return a (boost shared) pointer to it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_it Bucket::meshes_begin()
{
  return meshes_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_const_it Bucket::meshes_begin() const
{
  return meshes_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_it Bucket::meshes_end()
{
  return meshes_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the meshes_ map
//*******************************************************************|************************************************************//
Mesh_const_it Bucket::meshes_end() const
{
  return meshes_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_system(SystemBucket_ptr system, 
                                         const std::string &name)
{
  SystemBucket_it s_it = systems_.find(name);                        // check if a system with this name already exists
  if (s_it != systems_end())
  {
    dolfin::error(
            "SystemBucket named \"%s\" already exists in bucket",    // if it does, issue an error
                                  name.c_str());
  }
  else
  {
    systems_[name] = system;                                         // if not insert it into the systems_ map
    orderedsystems_[(int) systems_.size()] = system;                 // also insert it in the order it was registered in the 
                                                                     // orderedsystems_ map
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
SystemBucket_ptr Bucket::fetch_system(const std::string &name)
{
  SystemBucket_it s_it = systems_.find(name);                        // check if a system with this name exists in the bucket
  if (s_it == systems_end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "SystemBucket named \"%s\" does not exist in bucket.", 
                                name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return a pointer to it
  }
}

//*******************************************************************|************************************************************//
// return a constant (boost shared) pointer to a system bucket in the bucket data maps
//*******************************************************************|************************************************************//
const SystemBucket_ptr Bucket::fetch_system(const std::string &name) 
                                                              const
{
  SystemBucket_const_it s_it = systems_.find(name);                  // check if a system with this name exists in the bucket
  if (s_it == systems_end())
  {
    dolfin::error(
              "SystemBucket named \"%s\" does not exist in bucket.", // if it doesn't, throw an error
                                                      name.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return a pointer to it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_it Bucket::systems_begin()
{
  return systems_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_const_it Bucket::systems_begin() const
{
  return systems_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_it Bucket::systems_end()
{
  return systems_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the systems_ map
//*******************************************************************|************************************************************//
SystemBucket_const_it Bucket::systems_end() const
{
  return systems_.end();
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_it Bucket::orderedsystems_begin()
{
  return orderedsystems_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_const_it Bucket::orderedsystems_begin() const
{
  return orderedsystems_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_it Bucket::orderedsystems_end()
{
  return orderedsystems_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the orderedsystems_ map
//*******************************************************************|************************************************************//
int_SystemBucket_const_it Bucket::orderedsystems_end() const
{
  return orderedsystems_.end();
}

//*******************************************************************|************************************************************//
// register a base ufl symbol (associated with the derived ufl symbol) in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_baseuflsymbol(const std::string &baseuflsymbol, 
                                    const std::string &uflsymbol)
{
  string_it s_it = baseuflsymbols_.find(uflsymbol);                  // check if this ufl symbol already exists
  if (s_it != baseuflsymbols_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
            "Name with ufl symbol \"%s\" already exists in system.", 
                                              uflsymbol.c_str());
  }
  else
  {
    baseuflsymbols_[uflsymbol] = baseuflsymbol;                      // if it doesn't, assign the baseuflsymbol to the maps
  }
}

//*******************************************************************|************************************************************//
// return a string containing the base ufl symbol for a given ufl symbol from the bucket data maps
//*******************************************************************|************************************************************//
const std::string Bucket::fetch_baseuflsymbol(
                                const std::string &uflsymbol) const
{
  string_const_it s_it = baseuflsymbols_.find(uflsymbol);            // check if this ufl symbol exists
  if (s_it == baseuflsymbols_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
            "Name with uflsymbol \"%s\" does not exist in system.", 
                                                uflsymbol.c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return the string containing the base ufl symbol
  }
}

//*******************************************************************|************************************************************//
// check if a ufl symbol exists in the baseuflsymbols_ bucket data map
//*******************************************************************|************************************************************//
const bool Bucket::contains_baseuflsymbol(
                              const std::string &uflsymbol) const
{
  string_const_it s_it = baseuflsymbols_.find(uflsymbol);
  return s_it != baseuflsymbols_.end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a function with the given ufl symbol in the bucket data maps
//*******************************************************************|************************************************************//
void Bucket::register_uflsymbol(GenericFunction_ptr function, 
                                const std::string &uflsymbol)
{
  GenericFunction_it g_it = uflsymbols_.find(uflsymbol);             // check if the ufl symbol already exists
  if (g_it != uflsymbols_.end())
  {
    dolfin::error(                                                   // if it does, issue an error
    "GenericFunction with ufl symbol \"%s\" already exists in system.", 
                                  uflsymbol.c_str());
  }
  else
  {
    uflsymbols_[uflsymbol] = function;                               // if not, register the pointer in the maps
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to the function associated with the given ufl symbol
//*******************************************************************|************************************************************//
GenericFunction_ptr Bucket::fetch_uflsymbol(
                                const std::string &uflsymbol) const
{
  GenericFunction_const_it g_it = uflsymbols_.find(uflsymbol);       // check if the ufl symbol exists
  if (g_it == uflsymbols_.end())
  {
    dolfin::error(                                                   // if it doesn't, issue an error
    "GenericFunction with uflsymbol \"%s\" does not exist in system.", 
                                          uflsymbol.c_str());
  }
  else
  {
    return (*g_it).second;                                           // if it does, return a pointer to the associated function
  }
}

// Return a string describing the contents of the bucket
const std::string Bucket::str() const 
{
  std::stringstream s;
  int indent = 1;
  s << "Bucket " << name() << std::endl;
  s << uflsymbols_str(indent);
  s << coefficientspaces_str(indent);
  s << meshes_str(indent);
  s << systems_str(indent);
  return s.str();
}

// Return a string describing the contents of meshes_
const std::string Bucket::meshes_str(int indent) const
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
const std::string Bucket::systems_str(int indent) const
{
  std::stringstream s;
  for ( SystemBucket_const_it s_it = systems_begin(); s_it != systems_end(); s_it++ )
  {
    s << (*(*s_it).second).str(indent);
  }
  return s.str();
}

void Bucket::uflsymbols_fill_()
{
  for (SystemBucket_const_it s_it = systems_begin(); s_it != systems_end(); s_it++)
  {
    SystemBucket_ptr system = (*s_it).second;
    register_uflsymbol((*system).function(), (*system).uflsymbol());
    register_uflsymbol((*system).oldfunction(), (*system).uflsymbol()+"_n");
    register_uflsymbol((*system).iteratedfunction(), (*system).uflsymbol()+"_i");
    for (FunctionBucket_const_it f_it = (*system).fields_begin(); f_it != (*system).fields_end(); f_it++)
    {
      FunctionBucket_ptr field = (*f_it).second;
      register_uflsymbol((*system).function(), (*field).uflsymbol());
      register_uflsymbol((*system).oldfunction(), (*field).uflsymbol()+"_n");
      register_uflsymbol((*system).iteratedfunction(), (*field).uflsymbol()+"_i");
    }
    for (FunctionBucket_const_it f_it = (*system).coeffs_begin(); f_it != (*system).coeffs_end(); f_it++)
    {
      FunctionBucket_ptr coeff = (*f_it).second;
      register_uflsymbol((*coeff).function(), (*coeff).uflsymbol());
      register_uflsymbol((*coeff).oldfunction(), (*coeff).uflsymbol()+"_n");
      register_uflsymbol((*coeff).iteratedfunction(), (*coeff).uflsymbol()+"_i");
    }
  }
}

// Register a coefficientspace in the system
void Bucket::register_coefficientspace(FunctionSpace_ptr coefficientspace, std::string uflsymbol)
{
  // First check if a field with this name already exists
  FunctionSpace_it f_it = coefficientspaces_.find(uflsymbol);
  if (f_it != coefficientspaces_.end())
  {
    // if it does, issue an error
    dolfin::error("FunctionSpace with uflsymbol \"%s\" already exists in system coefficientspaces.", uflsymbol.c_str());
  }
  else
  {
    // if not then insert it into the maps
    coefficientspaces_[uflsymbol] = coefficientspace;
  }
}

// Check if the system contains a coefficientspace with the given uflsymbol
const bool Bucket::contains_coefficientspace(std::string uflsymbol) const
{
  // First check if a field with this uflsymbol already exists
  FunctionSpace_const_it f_it = coefficientspaces_.find(uflsymbol);
  return f_it != coefficientspaces_.end();
}

FunctionSpace_ptr Bucket::fetch_coefficientspace(std::string uflsymbol) const
{
  // First check if a functionspace with this uflsymbol already exists
  FunctionSpace_const_it f_it = coefficientspaces_.find(uflsymbol);
  if (f_it == coefficientspaces_.end())
  {
    // if it doesn't, issue an error
    std::cerr << coefficientspaces_str();
    dolfin::error("FunctionSpace with uflsymbol \"%s\" doesn't exist in system coefficientspaces.", uflsymbol.c_str());
  }
  else
  {
    // if it does return it
    return (*f_it).second;
  }
}

void Bucket::output()
{
  for (SystemBucket_it s_it = systems_begin(); s_it != systems_end(); s_it++)
  {
    (*(*s_it).second).output();
  } 
}

// Return a string describing the contents of uflsymbols_
const std::string Bucket::uflsymbols_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( GenericFunction_const_it g_it = uflsymbols_.begin(); g_it != uflsymbols_.end(); g_it++ )
  {
    if ((*g_it).second)
    {
      s << indentation << "UFLSymbol " << (*g_it).first << " associated" << std::endl;
    }
    else
    {
      s << indentation << "UFLSymbol " << (*g_it).first << " not associated" << std::endl;
    }
  }
  return s.str();
}

// Return a string describing the contents of coefficientspaces_
const std::string Bucket::coefficientspaces_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( FunctionSpace_const_it f_it = coefficientspaces_.begin(); f_it != coefficientspaces_.end(); f_it++ )
  {
    s << indentation << "CoefficientSpace for " << (*f_it).first  << std::endl;
  }
  return s.str();
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
//    dolfin::error("SystemBucket named \"%s\" already exists in bucket.", name.c_str());
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

