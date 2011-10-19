
#ifndef __PYTHON_INSTANCE_H
#define __PYTHON_INSTANCE_H

#include "Python.h"
#include <string>

namespace buckettools
{
  //*****************************************************************|************************************************************//
  // PythonInstance class:
  //
  // The PythonInstance class wraps python functionality so that it is accessible from C++
  // NOTE: This assumes a function name val is available and runs that and only that!
  //*****************************************************************|************************************************************//
  class PythonInstance
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    PythonInstance(const std::string &function);                     // specific constructor (takes a string with the python function)
    
    ~PythonInstance();                                               // default destructor
    
    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string function() const                               // return a constant string containing the function
    { return function_; }
    
    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    PyObject* call(PyObject *pArgs) const;                           // run the function contained in this python instance

    const int number_arguments() const                               // return the number of arguments expected by this pythoninstance
    { return nargs_; }

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    const std::string function_;                                     // the python function string

    PyObject *pMain_, *pGlobals_, *pLocals_, *pCode_, *pFunc_;       // python objects used to run the function (and cacheable between calls)

    int nargs_;                                                      // the number of arguments this python function takes
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_();                                                    // initialize the python instance
    
    //***************************************************************|***********************************************************//
    // Clean up
    //***************************************************************|***********************************************************//

    void clean_();                                                   // clean the python instance
    
  };
}
#endif
