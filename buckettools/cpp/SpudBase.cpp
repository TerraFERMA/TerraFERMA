#include "BoostTypes.h"
#include "SpudBase.h"
#include "PythonExpression.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

//*******************************************************************|************************************************************//
// initialize an expression from an optionpath assuming the buckettools schema
//*******************************************************************|************************************************************//
Expression_ptr buckettools::initialize_expression(const std::string 
                                                        &optionpath, 
                                                  const int &size, 
                                                  const std::vector<int> 
                                                        &shape)
{
  Spud::OptionError serr;                                            // spud option error
  Expression_ptr expression;                                         // declare the pointer that will be returned
  
  std::stringstream constbuffer, pybuffer;                           // some string assembly streams
  constbuffer.str(""); constbuffer << optionpath << "/constant";     // for a constant
  pybuffer.str(""); pybuffer << optionpath << "/python";             // for a python function
  
  if (Spud::have_option(constbuffer.str()))                          // constant requested?
  {
    int rank;
    serr = Spud::get_option_rank(constbuffer.str(), rank);           // find out the rank
    spud_err(constbuffer.str(), serr);
    if(rank==0)                                                      // scalar
    {
      double value;
      serr = Spud::get_option(constbuffer.str(), value); 
      spud_err(constbuffer.str(), serr);
      expression.reset(new dolfin::Constant(value));
    }
    else if (rank==1)                                                // vector
    {
      std::vector<double> values;
      serr = Spud::get_option(constbuffer.str(), values); 
      spud_err(constbuffer.str(), serr);
      assert(values.size()==size);
      expression.reset(new dolfin::Constant(values));
    }
//    else if (rank==2)
//    {
//      std::vector<int> value_shape;
//      std::vector<double> values; // not sure this will work 
//                                  // (might have to be std::vector< std::vector<double> >
//                                  //  but this disagrees with the DOLFIN interface)
//      serr = Spud::get_option_shape(constbuffer.str(), value_shape); spud_err(constbuffer.str(), serr);
//      serr = Spud::get_option(constbuffer.str(), values); spud_err(constbuffer.str(), serr);
//      expression.reset(new dolfin::Constant());
//    }
//    else
//    {
//      dolfin::error("Unknown rank in init_exp_");
//    }
    else                                                             // unable to deal with this rank
    {
      dolfin::error("Don't deal with rank > 1 yet.");
    }
  } 
  else if (Spud::have_option(pybuffer.str()))                        // python requested
  {
    std::string pyfunction;                                          // the python function string
    serr = Spud::get_option(pybuffer.str(), pyfunction); 
    spud_err(pybuffer.str(), serr);
    
                                                                     // rank of a python function isn't in the default spud base
                                                                     // language so have added it... but it comes out as a string 
                                                                     // of course!
    std::stringstream buffer;
    std::string rankstring;                                          // bit of a hack
    buffer.str(""); buffer << pybuffer.str() << "/rank";
    serr = Spud::get_option(buffer.str(), rankstring); 
    spud_err(buffer.str(), serr);
    
    int rank;
    rank = atoi(rankstring.c_str());
    if(rank==0)                                                      // scalar
    {
      expression.reset(new buckettools::PythonExpression(pyfunction));
    }
    else if (rank==1)                                                // vector
    {
      expression.reset(new buckettools::PythonExpression(size, pyfunction));
    }
    else
    {
      dolfin::error("Don't deal with rank > 1 yet.");
    }
  }
  else                                                               // unknown expression type
  {
    dolfin::error("Unknown way of specifying expression.");
  }
  
  return expression;                                                 // return the initialized expression
  
}


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
