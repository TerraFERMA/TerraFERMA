#ifndef __VISUALIZATION_WRAPPER_H
#define __VISUALIZATION_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  FunctionSpace_ptr ufc_fetch_visualization_functionspace(           // return a (boost shared) pointer to a functionspace for
                                      const std::string &meshname,   // visualization given a (boost shared) pointer to a mesh and
                                      Mesh_ptr mesh);                // its name

}

#endif
