
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  
  class System
  {
  private:

    std::string name_;

  public:
    
    System(Mesh_ptr mesh)
    { System("uninitialised_name", mesh); }

    System(std::string name, Mesh_ptr mesh);
    
    ~System();

    std::string name()
    { return name_; }

  protected:
    
    Mesh_ptr mesh_;

    FunctionSpace_ptr functionspace_;
    
  };

  typedef boost::shared_ptr< System > System_ptr;

}
#endif
