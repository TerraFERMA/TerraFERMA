
#include "DolfinBoostTypes.h"
#include "System.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

System::System(std::string name, Mesh_ptr mesh) : name_(name), mesh_(mesh)
{
  // Do nothing
}

System::~System()
{
  // Do nothing
}

// Register a dolfin subfunctionspace
void System::register_subfunctionspace(FunctionSpace_ptr subfunctionspace, std::string name)
{
  // First check if a mesh with this name already exists
  FunctionSpace_it f_it = subfunctionspaces_.find(name);
  if (f_it != subfunctionspaces_.end())
  {
    // if it does, issue an error
    dolfin::error("Subfunctionspace named \"%s\" already exists in system.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    subfunctionspaces_[name] = subfunctionspace;
  }
}

