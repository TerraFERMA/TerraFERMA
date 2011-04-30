
#ifndef __PYTHON_INSTANCE_H
#define __PYTHON_INSTANCE_H

#include "Python.h"
#include <string>

namespace buckettools
{
  class PythonInstance
  {
  private:
    const std::string function_;
    PyObject *pMain_, *pGlobals_, *pLocals_, *pCode_, *pFunc_;
    
    void init_();
    
    void clean_();
    
  public:
    PythonInstance(std::string function);
    
    ~PythonInstance();
    
    std::string function() const;
    
    PyObject* call(PyObject *pArgs) const;
  };
}

#endif