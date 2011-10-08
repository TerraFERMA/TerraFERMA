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
                           const Spud::OptionError &error)
{
  switch (error)
  {
    case Spud::SPUD_NO_ERROR:                                        // no error
      break;
    case Spud::SPUD_KEY_ERROR:                                       // key error, e.g. unknown key
      dolfin::error("Option key error. Key is: %s", 
                                                optionpath.c_str());
      break;
    case Spud::SPUD_TYPE_ERROR:                                      // type error, e.g. integer not float
      dolfin::error("Option type error. Key is: %s", 
                                                optionpath.c_str());
      break;
    case Spud::SPUD_RANK_ERROR:                                      // rank error, e.g. wrong rank passed in
      dolfin::error("Option rank error. Key is: %s", 
                                                optionpath.c_str());
      break;
    case Spud::SPUD_SHAPE_ERROR:                                     // shape error, e.g. wrong shape
      dolfin::error("Option shape error. Key is: %s", 
                                                optionpath.c_str());
      break;
    case Spud::SPUD_FILE_ERROR:                                      // file error
      dolfin::error("Option file error. Filename is: %s", 
                                                optionpath.c_str());
      break;
    case Spud::SPUD_NEW_KEY_WARNING:                                 // adding a new key
      dolfin::error(                                                 // FIXME: not an error
                "Option warning. Key is not in the options tree: %s",
                                                optionpath.c_str());
      break;
    case Spud::SPUD_ATTR_SET_FAILED_WARNING:                         // attribute error
      dolfin::error(                                                 // FIXME: not an error
      "Option warning. Option cannot be set as an attribute. Key is %s",
                                                optionpath.c_str());
      break;
    default:                                                         // unknown error
      dolfin::error("Unknown option error. Key is: ", 
                                                optionpath.c_str());
  }

}
