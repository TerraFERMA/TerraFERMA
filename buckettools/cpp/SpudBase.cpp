#include "BoostTypes.h"
#include "SpudBase.h"
#include "PythonExpression.h"
#include <dolfin.h>
#include <string>
#include <spud>

using namespace buckettools;

// initialize an expression using spud
Expression_ptr buckettools::initialize_expression(const std::string &optionpath, const int &size, const std::vector<int> &shape)
{
  Spud::OptionError serr;
  // declare the pointer we'll be returning
  Expression_ptr expression;
  
  std::stringstream constbuffer, pybuffer;
  constbuffer.str(""); constbuffer << optionpath << "/constant";
  pybuffer.str(""); pybuffer << optionpath << "/python";
  
  // Are we constant or python (or something else we don't know about yet)?
  if (Spud::have_option(constbuffer.str()))
  {
    // constant
    int rank;
    serr = Spud::get_option_rank(constbuffer.str(), rank); spud_err(constbuffer.str(), serr);
    if(rank==0)
    {
      double value;
      serr = Spud::get_option(constbuffer.str(), value); spud_err(constbuffer.str(), serr);
      expression.reset(new dolfin::Constant(value));
    }
    else if (rank==1)
    {
      std::vector<double> values;
      serr = Spud::get_option(constbuffer.str(), values); spud_err(constbuffer.str(), serr);
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
    else
    {
      dolfin::error("Don't deal with rank > 1 yet.");
    }
  } 
  else if (Spud::have_option(pybuffer.str()))
  {
    // python
    std::string pyfunction;
    serr = Spud::get_option(pybuffer.str(), pyfunction); spud_err(pybuffer.str(), serr);
    
    // rank of a python function isn't in the default spud base language
    // so have added it... but it comes out as a string of course!
    std::stringstream buffer;
    std::string rankstring; // bit of a hack
    buffer.str(""); buffer << pybuffer.str() << "/rank";
    serr = Spud::get_option(buffer.str(), rankstring); spud_err(buffer.str(), serr);
    
    int rank;
    rank = atoi(rankstring.c_str());
    if(rank==0)
    {
      expression.reset(new buckettools::PythonExpression(pyfunction));
    }
    else if (rank==1)
    {
      expression.reset(new buckettools::PythonExpression(size, pyfunction));
    }
    else
    {
      dolfin::error("Don't deal with rank > 1 yet.");
    }
  }
  else
  {
    dolfin::error("Unknown way of specifying expression.");
  }
  
  return expression;
  
}


// Handle option errors through dolfin::error
void buckettools::spud_err(const std::string &optionpath, 
                           const Spud::OptionError &error)
{
  switch (error)
  {
    case Spud::SPUD_NO_ERROR:
      break;
    case Spud::SPUD_KEY_ERROR:
      dolfin::error("Option key error. Key is: %s", optionpath.c_str());
      break;
    case Spud::SPUD_TYPE_ERROR:
      dolfin::error("Option type error. Key is: ", optionpath.c_str());
      break;
    case Spud::SPUD_RANK_ERROR:
      dolfin::error("Option rank error. Key is: ", optionpath.c_str());
      break;
    case Spud::SPUD_SHAPE_ERROR:
      dolfin::error("Option shape error. Key is: ", optionpath.c_str());
      break;
    case Spud::SPUD_FILE_ERROR:
      dolfin::error("Option file error. Filename is: ", optionpath.c_str());
      break;
    case Spud::SPUD_NEW_KEY_WARNING:
      dolfin::error("Option warning. Key is not in the options tree: ", optionpath.c_str());
      break;
    case Spud::SPUD_ATTR_SET_FAILED_WARNING:
      dolfin::error("Option warning. Option cannot be set as an attribute. Key is ", optionpath.c_str());
      break;
    default:
      dolfin::error("Unknown option error. Key is: ", optionpath.c_str());
  }

}
