
#ifndef __SYSTEMBUCKET_H
#define __SYSTEMBUCKET_H

#include "Detectors.h"
#include "PythonDetectors.h"
#include <dolfin.h>

namespace buckettools
{
  
  class SystemBucket
  {
  private:
    
    std::map< std::string, FunctionSpace_ptr > functionspaces_;
    
    std::map< std::string, DirichletBC_ptr > dirichletbcs_;
    
    std::vector< GenericFunction_ptr > bcexps_;
    
    std::map< std::string, Function_ptr > functions_;
    
    std::map< std::string, Detectors_ptr > detectors_;
    
    
    void clean_();
    
  public:
    
    SystemBucket();
    
    ~SystemBucket();
    
    void register_mesh(Mesh_ptr mesh, std::string name);
    
    void register_meshfunction(MeshFunction_uint_ptr meshfunction, std::string name);
    
    void register_functionspace(FunctionSpace_ptr functionspace, std::string name);
    
    void register_dirichletbc(DirichletBC_ptr dirichletbc, std::string name);
    
    void register_detector(Detectors_ptr detector, std::string name);
    
    void register_function(Function_ptr function, std::string name);
    
    Mesh_ptr fetch_mesh(const std::string name);
    
    MeshFunction_uint_ptr fetch_meshfunction(const std::string name);
    
    FunctionSpace_ptr fetch_functionspace(const std::string name);
    
//     Function_ptr fetch_function(const std::string name);
    
    Function_ptr fetch_function(const std::string name);
    
    Detectors_ptr fetch_detector(const std::string name);
    
    std::map< std::string, Detectors_ptr >::iterator detectors_begin();
    
    std::map< std::string, Detectors_ptr >::const_iterator detectors_begin() const;
    
    std::map< std::string, Detectors_ptr >::iterator detectors_end();
    
    std::map< std::string, Detectors_ptr >::const_iterator detectors_end() const;
    
    std::map< std::string, DirichletBC_ptr >::iterator dirichletbcs_begin();
    
    std::map< std::string, DirichletBC_ptr >::const_iterator dirichletbcs_begin() const;
    
    std::map< std::string, DirichletBC_ptr >::iterator dirichletbcs_end();
    
    std::map< std::string, DirichletBC_ptr >::const_iterator dirichletbcs_end() const;
    
    std::map< std::string, Function_ptr >::iterator functions_begin();
    
    std::map< std::string, Function_ptr >::const_iterator functions_begin() const;
    
    std::map< std::string, Function_ptr >::iterator functions_end();
    
    std::map< std::string, Function_ptr >::const_iterator functions_end() const;
    
  protected:
    
    void register_bcexp(GenericFunction_ptr bcexp);
    
  };
}
#endif
