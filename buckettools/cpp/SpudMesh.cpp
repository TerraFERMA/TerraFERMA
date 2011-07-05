#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor for spudbucket derived class
SpudMesh::SpudMesh(std::string name, std::string option_path, std::string filename) : name_(name), option_path_(option_path), dolfin::Mesh(filename)
{
  // Do nothing
}

// Default destructor for spudbucket derived class
SpudMesh::~SpudMesh()
{
  // Do nothing
}

