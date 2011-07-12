
#ifndef __PYTHON_INSTANCE_H
#define __PYTHON_INSTANCE_H

#include "Python.h"
#include <string>

namespace buckettools
{
  // The PythonInstance class wraps python functionality so that it is
  // accessible from C++
  // NOTE: This assumes a function name val is available and runs that
  // and only that!
  class PythonInstance
  {
  // only accessible to this class
  private:

    // a string describing a python function
    const std::string function_;

    // Python objects for the function
    PyObject *pMain_, *pGlobals_, *pLocals_, *pCode_, *pFunc_;
    
    // Initialize the python instance
    void init_();
    
    // Clean up the python instance
    void clean_();
    
  // accessible to everyone
  public:

    // No default constructor as we require a python function

    // Specific constructor
    PythonInstance(const std::string &function);
    
    // Default destructor
    ~PythonInstance();
    
    // Return the function
    std::string function() const
    { return function_; }
    
    // Run the python instance by calling the function
    PyObject* call(PyObject *pArgs) const;
  };
}

#endif
