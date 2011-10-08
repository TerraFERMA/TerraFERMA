
#ifndef __SPUD_BASE_H
#define __SPUD_BASE_H

#include "BoostTypes.h"
#include <dolfin.h>
#include <spud>

namespace buckettools
{
  //*****************************************************************|************************************************************//
  // A collection of tools that are spud (and schema) specific but do not depend on any class
  // objects in the bucketools library
  //*****************************************************************|************************************************************//

  void spud_err(const std::string &optionpath,                       // translate the spud error codes into dolfin errors
                const Spud::OptionError &error);

}

#endif
