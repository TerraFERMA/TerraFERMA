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

#include "BoostTypes.h"
#include "SpudBase.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// return dolfin errors after failures in spud
//*******************************************************************|************************************************************//
void buckettools::spud_err(const std::string &optionpath, 
                           const Spud::OptionError &error,
                           const Spud::OptionError accept)
{
  switch (error)
  {
    case Spud::SPUD_NO_ERROR:                                        // no error
      break;
    case Spud::SPUD_KEY_ERROR:                                       // key error, e.g. unknown key
      if(accept!=Spud::SPUD_KEY_ERROR)
      {
        dolfin::error("Option key error. Key is: %s", 
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_TYPE_ERROR:                                      // type error, e.g. integer not float
      if(accept!=Spud::SPUD_TYPE_ERROR)
      {
        dolfin::error("Option type error. Key is: %s", 
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_RANK_ERROR:                                      // rank error, e.g. wrong rank passed in
      if(accept!=Spud::SPUD_RANK_ERROR)
      {
        dolfin::error("Option rank error. Key is: %s", 
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_SHAPE_ERROR:                                     // shape error, e.g. wrong shape
      if(accept!=Spud::SPUD_SHAPE_ERROR)
      {
        dolfin::error("Option shape error. Key is: %s", 
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_FILE_ERROR:                                      // file error
      if(accept!=Spud::SPUD_FILE_ERROR)
      {
        dolfin::error("Option file error. Filename is: %s", 
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_NEW_KEY_WARNING:                                 // adding a new key
      if(accept!=Spud::SPUD_NEW_KEY_WARNING)
      {
        dolfin::warning(
                  "Option warning. Key is not in the options tree: %s",
                                                 optionpath.c_str());
      }
      break;
    case Spud::SPUD_ATTR_SET_FAILED_WARNING:                         // attribute error
      if(accept!=Spud::SPUD_ATTR_SET_FAILED_WARNING)
      {
        dolfin::warning(
        "Option warning. Option cannot be set as an attribute. Key is %s",
                                                 optionpath.c_str());
      }
      break;
    default:                                                         // unknown error
      dolfin::error("Unknown option error. Key is: ", 
                                                 optionpath.c_str());
  }

}
