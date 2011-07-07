
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
