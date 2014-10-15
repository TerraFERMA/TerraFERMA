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
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// return TF errors after failures in spud
//*******************************************************************|************************************************************//
void buckettools::spud_error(const std::string &optionpath, 
                             const Spud::OptionError &serr,
                             const std::string &filename, 
                             const int &line,
                             const Spud::OptionError &accept)
{
  switch (serr)
  {
    case Spud::SPUD_NO_ERROR:                                        // no error
      break;
    case Spud::SPUD_KEY_ERROR:                                       // key error, e.g. unknown key
      if(accept!=Spud::SPUD_KEY_ERROR)
      {
        error(filename, line, "Spud option key error.", "Key is: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_TYPE_ERROR:                                      // type error, e.g. integer not float
      if(accept!=Spud::SPUD_TYPE_ERROR)
      {
        error(filename, line, "Spud option key error.", "Key is: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_RANK_ERROR:                                      // rank error, e.g. wrong rank passed in
      if(accept!=Spud::SPUD_RANK_ERROR)
      {
        error(filename, line, "Spud option rank error.", "Key is: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_SHAPE_ERROR:                                     // shape error, e.g. wrong shape
      if(accept!=Spud::SPUD_SHAPE_ERROR)
      {
        error(filename, line, "Spud option shape error.", "Key is: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_FILE_ERROR:                                      // file error
      if(accept!=Spud::SPUD_FILE_ERROR)
      {
        error(filename, line, "Spud option file error.", "Key is: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_NEW_KEY_WARNING:                                 // adding a new key
      if(accept!=Spud::SPUD_NEW_KEY_WARNING)
      {
        warning(filename, line, "Spud option warning.", "Key is not in the options tree: %s", optionpath.c_str());
      }
      break;
    case Spud::SPUD_ATTR_SET_FAILED_WARNING:                         // attribute error
      if(accept!=Spud::SPUD_ATTR_SET_FAILED_WARNING)
      {
        warning(filename, line, "Spud option warning.", "Option cannot be set as an attribute. Key is: %s", optionpath.c_str());
      }
      break;
    default:                                                         // unknown error
      error(filename, line, "Spud unknown option error.", "Spud error code is: %d. Key is: %s", serr, optionpath.c_str());
  }

}
