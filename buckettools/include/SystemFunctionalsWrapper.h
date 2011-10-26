#ifndef __SYSTEMFUNCTIONALS_WRAPPER_H
#define __SYSTEMFUNCTIONALS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(
                                      const std::string              // return a (boost shared) pointer to a functionspace for a 
                                            &systemname,             // coefficient function from a functional given a mesh, a
                                      const std::string              // function name, a functional name and a (base) ufl symbol
                                            &functionname, 
                                      const std::string 
                                            &funcionalname, 
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh);

  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(
                                      const std::string              // return a (boost shared) pointer to a functionspace for a 
                                            &systemname,             // coefficient function from a constant functional given a mesh, a
                                      const std::string              // function name and a (base) ufl symbol
                                            &coefficientname, 
                                      const std::string &uflsymbol, 
                                      Mesh_ptr mesh);

  Form_ptr ufc_fetch_functional(const std::string &systemname,       // return a (boost shared) pointer to a form from a functional
                          const std::string &functionname,           // given a mesh, a function name and a functional name
                          const std::string &functionalname, 
                          Mesh_ptr mesh);

  Form_ptr ufc_fetch_functional(const std::string &systemname,       // return a (boost shared) pointer to a form for a constant
                          const std::string &functionname,           // from a functional given a mesh and a function name
                          Mesh_ptr mesh);

}

#endif
