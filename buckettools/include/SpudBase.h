
#ifndef __SPUD_BASE_H
#define __SPUD_BASE_H

#include "DolfinBoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  
  Expression_ptr initialize_expression(const std::string &optionpath, const uint &dimension);

}

#endif
