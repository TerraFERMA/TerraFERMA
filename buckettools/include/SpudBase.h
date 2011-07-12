
#ifndef __SPUD_BASE_H
#define __SPUD_BASE_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  // A collection of tools that are spud (and schema) specific but do not depend on any class
  // objects in the bucketools library

  // Initialize an expression based on a optionpath (and dimension for python vectors)
  Expression_ptr initialize_expression(const std::string &optionpath, const uint &dimension);

}

#endif
