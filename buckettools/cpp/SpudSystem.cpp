
#include "DolfinBoostTypes.h"
#include "SpudSystem.h"
#include "System.h"
#include "SystemsWrapper.h"
#include <dolfin.h>
#include <string>
#include <spud.h>

using namespace buckettools;

SpudSystem::SpudSystem(std::string name, std::string optionpath, Mesh_ptr mesh) : optionpath_(optionpath), System(name, mesh)
{
  // Do nothing
}

SpudSystem::~SpudSystem()
{
  // Do nothing
}

//void SpudSystem::fill(Bucket_ptr bucket)
void SpudSystem::fill()
{
  // A buffer to put option paths in
  std::stringstream buffer;

}

