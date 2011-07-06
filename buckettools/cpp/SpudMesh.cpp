#include <dolfin.h>
#include <string>
#include "SpudMesh.h"

using namespace buckettools;

// Default constructor for spudbucket derived class
SpudMesh::SpudMesh(std::string name, std::string optionpath, std::string filename) : name_(name), optionpath_(optionpath), dolfin::Mesh(filename)
{
  // Do nothing
}

// Default destructor for spudbucket derived class
SpudMesh::~SpudMesh()
{
  // Do nothing
}

