
#ifndef __SPUD_BUCKET_H
#define __SPUD_BUCKET_H

#include "Bucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  class SpudBucket : public Bucket
  {
  private:
    
    std::string optionpath_;

    std::map< std::string, std::string > mesh_optionpaths_;

    void meshes_fill_(const std::string &optionpath);
    
    //void system_fill_(const std::string &optionpath);
    //
    //void detectors_fill_();
    //
    //void bc_fill_(const std::string bcpath, 
    //              const int funci,
    //              const int dimension,
    //              FunctionSpace_ptr subsysspace,
    //              const std::string sysname,
    //              const std::string funcname);
    //
    //GenericFunction_ptr init_exp_(const std::string path, const int dimension);
    
  public:
    
    SpudBucket()
    { SpudBucket("uninitialised_name", ""); }
    
    SpudBucket(std::string name)
    { SpudBucket(name, ""); }
    
    SpudBucket(std::string name, std::string option_path);
    
    void register_mesh(Mesh_ptr mesh, std::string name, std::string optionpath);

    std::string fetch_mesh_optionpath(const std::string name);

    std::string meshes_str() const;

    virtual ~SpudBucket();
    
    void fill();
    
    string_it mesh_optionpaths_begin();

    string_const_it mesh_optionpaths_begin() const;

    string_it mesh_optionpaths_end();

    string_const_it mesh_optionpaths_end() const;
 
  };
}
#endif
