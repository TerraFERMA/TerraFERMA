
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

    Mesh_ptr mesh_;

    FunctionSpace_ptr functionspace_;

  public:
    
    System(Mesh_ptr mesh)
    { System("uninitialised_name", mesh); }

    System(std::string name, Mesh_ptr mesh);
    
    ~System();
    
  protected:
    
    
  };

  typedef boost::shared_ptr< System > System_ptr;

}
#endif
