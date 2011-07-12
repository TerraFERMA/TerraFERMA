#include "DolfinBoostTypes.h"
#include "SpudBase.h"
#include <dolfin.h>
#include <string>
#include <spud.h>

using namespace buckettools;

Expression_ptr initialize_expression(const std::string &optionpath, const uint &dimension)
{
  Expression_ptr expression;
  
  std::stringstream constbuffer, pybuffer;
  constbuffer.str(""); constbuffer << optionpath << "/constant";
  pybuffer.str(""); pybuffer << optionpath << "/python";
  
  if (Spud::have_option(constbuffer.str()))
  {
    int rank;
    Spud::get_option_rank(constbuffer.str(), rank);
    if(rank==0)
    {
      double value;
      Spud::get_option(constbuffer.str(), value);
      expression.reset(new dolfin::Constant(value));
    }
    else if (rank==1)
    {
      std::vector<double> values;
      Spud::get_option(constbuffer.str(), values);
      expression.reset(new dolfin::Constant(values));
    }
//    else if (rank==2)
//    {
//      std::vector<int> value_shape;
//      std::vector<double> values; // not sure this will work 
//                                  // (might have to be std::vector< std::vector<double> >
//                                  //  but this disagrees with the DOLFIN interface)
//      Spud::get_option_shape(constbuffer.str(), value_shape);
//      Spud::get_option(constbuffer.str(), values);
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
    std::string pyfunction;
    Spud::get_option(pybuffer.str(), pyfunction);
    
    // rank of a python function isn't in the default spud base language
    // so have added it... but it comes out as a string of course!
    std::stringstream buffer;
    std::string rankstring; // bit of a hack
    buffer.str(""); buffer << pybuffer.str() << "/rank";
    Spud::get_option(buffer.str(), rankstring);
    
    int rank;
    rank = atoi(rankstring.c_str());
    if(rank==0)
    {
      expression.reset(new buckettools::PythonExpression(pyfunction));
    }
    else if (rank==1)
    {
      expression.reset(new buckettools::PythonExpression(dimension, pyfunction));
    }
    else
    {
      dolfin::error("Don't deal with rank > 1 yet.");
    }
  }
  else
  {
    dolfin::error("Unknown way of specifying bc expression.");
  }
  
  return expression;
  
}
