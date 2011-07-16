#ifndef __SYSTEMS_WRAPPER_H
#define __SYSTEMS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{
  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).
  FunctionSpace_ptr ufc_fetch_functionspace(std::string systemname, Mesh_ptr mesh);

  // A function to return a functionspace from a system given a mesh and a solvername.
  FunctionSpace_ptr ufc_fetch_functionspace(std::string systemname, std::string solvername, Mesh_ptr mesh);

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.
  FunctionSpace_ptr ufc_fetch_coefficientspace(std::string systemname, std::string solvername, std::string coefficientname, Mesh_ptr mesh);

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.
  FunctionSpace_ptr ufc_fetch_coefficientspace(std::string systemname, std::string functionname, std::string funcionalname, std::string coefficientname, Mesh_ptr mesh);

  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.
  Form_ptr ufc_fetch_form(std::string systemname, std::string solvername, std::string solvertype, std::string formname, FunctionSpace_ptr functionspace);

  // A function to return a functional from a system-function set given a mesh and a functionalname.
  Form_ptr ufc_fetch_functional(std::string systemname, std::string functionname, std::string functionalname, Mesh_ptr mesh);

}

#endif
