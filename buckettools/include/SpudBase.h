
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

  Expression_ptr initialize_expression(const std::string &optionpath,// initialize an expression based on an optionpath (and a size
                                       const int *size,              // and shape for vectors and tensors)
                                       const std::vector<int>        // NOTE: this declaration has to come first!
                                                            *shape);

  Expression_ptr initialize_expression(const std::string &optionpath);// initialize an expression based on an optionpath

  Expression_ptr initialize_expression(const std::string &optionpath,// initialize an expression based on an optionpath (and a size
                                       const int *size);              // and shape for vectors)

  Expression_ptr initialize_expression(const std::string &optionpath,// initialize an expression based on an optionpath (and a
                                       const std::vector<int> *shape);// and shape for tensors)

  void spud_err(const std::string &optionpath,                       // translate the spud error codes into dolfin errors
                const Spud::OptionError &error);

}

#endif
