#ifndef __SYSTEMEXPRESSIONS_WRAPPER_H
#define __SYSTEMEXPRESSIONS_WRAPPER_H

#include "BoostTypes.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // Header definitions for code that is automatically generated on a per options file basis.
  // These functions provide an interface between the bucket and the ufc, which is model specific.
  // DO NOT CHANGE THESE INTERFACES WITHOUT UPDATING THE CODE GENERATION SCRIPT.
  //*****************************************************************|************************************************************//

  Expression_ptr cpp_fetch_expression(const std::string &systemname, // return a (boost shared) pointer to an expression 
                          const std::string &functionname,           // given a system, function & expression name and a bucket
                          const std::string &expressiontype,         // and an expression type (initial_condition, boudary_condition,
                          const std::string &expressionname,         // value)
                          const uint &size, 
                          const std::vector<int> &shape, 
                          const Bucket *bucket,  
                          const SystemBucket *system,
                          const double_ptr time);

  void cpp_init_expression(Expression_ptr expression,                // initialize the given (boost shared) pointer to an expression
                          const std::string &systemname,             // given a system, function & expression name & an expression type
                          const std::string &functionname,
                          const std::string &expressiontype,
                          const std::string &expressionname);

}

#endif
