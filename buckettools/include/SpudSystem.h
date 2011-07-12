
#ifndef __SPUD_SYSTEM_H
#define __SPUD_SYSTEM_H

#include "DolfinBoostTypes.h"
#include "System.h"
#include <dolfin.h>

namespace buckettools
{
  
  class SpudSystem : public System
  {
  private:

    std::string optionpath_;

    void fields_fill_(const std::string &optionpath, const int &field_i);

    void bc_fill_(const std::string &optionpath, 
                  const int &field_i,
                  const int dimension,
                  FunctionSpace_ptr subsysspace,
                  const std::string &fieldname);
    
  public:
    
    //SpudSystem(Mesh_ptr mesh)
    //{ SpudSystem("uninitialised_name", "uninitialised_path", mesh); }

    //SpudSystem(std::string name, Mesh_ptr mesh)
    //{ SpudSystem(name, "uninitialised_path", mesh); }

    SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh);
    
    ~SpudSystem();

    void fill();

    std::string optionpath()
    { return optionpath_; }
    
  protected:
    
    
  };
 
  typedef boost::shared_ptr< SpudSystem > SpudSystem_ptr;

}
#endif
