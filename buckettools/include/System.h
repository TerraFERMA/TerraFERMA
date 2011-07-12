
#ifndef __SYSTEM_H
#define __SYSTEM_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{

  typedef std::map< std::string, FunctionSpace_ptr >::iterator        FunctionSpace_it;
  typedef std::map< std::string, FunctionSpace_ptr >::const_iterator  FunctionSpace_const_it;
  
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

    void register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name);

  protected:
    
    Mesh_ptr mesh_;

    FunctionSpace_ptr functionspace_;

    Function_ptr function_, oldfunction_, iteratedfunction_;

    std::map< std::string, FunctionSpace_ptr > subfunctionspaces_;
    
  };

  typedef boost::shared_ptr< System > System_ptr;

}
#endif
