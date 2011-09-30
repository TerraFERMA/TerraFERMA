#ifndef __SYSTEMS_WRAPPER_H
#define __SYSTEMS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{
  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  FunctionSpace_ptr ufc_fetch_functionspace(const std::string        // return a (boost shared) pointer to a functionspace from a
                                            &systemname,             // system given a (boost shared) pointer to a mesh (defaults to the
                                      Mesh_ptr mesh);                // first solver in the system as they should all be the same)

  FunctionSpace_ptr ufc_fetch_functionspace(const std::string        // return a (boost shared) pointer to a functionspace from a
                                            &systemname,             // system given a mesh and a solver name (i.e. does not default
                                      const std::string &solvername, // as above) 
                                      Mesh_ptr mesh);

  FunctionSpace_ptr ufc_fetch_coefficientspace(const std::string     // return a (boost shared) pointer to a functionspace for a
                                            &systemname,             // coefficient function from a solver given a mesh, a solver name
                                      const std::string &solvername, // and a (base) ufl symbol
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh);

  FunctionSpace_ptr ufc_fetch_coefficientspace(const std::string     // return a (boost shared) pointer to a functionspace for a 
                                            &systemname,             // coefficient function from a functional given a mesh, a
                                      const std::string              // function name, a functional name and a (base) ufl symbol
                                            &functionname, 
                                      const std::string 
                                            &funcionalname, 
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh);

  Form_ptr ufc_fetch_form(const std::string &systemname,             // return a (boost shared) pointer to a form from a solver
                          const std::string &solvername,             // given a functionspace, a solver name, a solver type
                          const std::string &solvertype,             // and a form name
                          const std::string &formname, 
                          FunctionSpace_ptr functionspace);

  Form_ptr ufc_fetch_functional(const std::string &systemname,       // return a (boost shared) pointer to a form from a functional
                          const std::string &functionname,           // given a mesh, a function name and a functional name
                          const std::string &functionalname, 
                          Mesh_ptr mesh);

  Form_ptr ufc_fetch_functional(const std::string &systemname,       // return a (boost shared) pointer to a form for a constant
                          const std::string &functionname,           // from a functional given a mesh and a function name
                          Mesh_ptr mesh);

}

#endif
