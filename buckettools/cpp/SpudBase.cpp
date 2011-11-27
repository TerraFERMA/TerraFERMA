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
