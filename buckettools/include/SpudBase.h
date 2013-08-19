// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


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
                const Spud::OptionError &error,
                const Spud::OptionError accept=Spud::SPUD_NO_ERROR);

}

#endif
