
#include "PythonInstance.h"
#include "Python.h"
#include <string>
#include <dolfin.h>

using namespace buckettools;

// Specific constructor
PythonInstance::PythonInstance(const std::string &function) : function_(function)
{
  init_();
}

// Default destructor
PythonInstance::~PythonInstance()
{
  clean_();
}

// Initialize the python objects contained in the python instance
void PythonInstance::init_()
{
  if(not Py_IsInitialized())
  {
    Py_Initialize();
  }
  
  pMain_ = PyImport_AddModule("__main__");
  pGlobals_ = PyModule_GetDict(pMain_);
  pLocals_ = PyDict_New();
  
  pCode_ = PyRun_String((char*)function_.c_str(), Py_file_input, pGlobals_, pLocals_);
  
  pFunc_ = PyDict_GetItemString(pLocals_, "val");
  
  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    dolfin::error("In PythonInstance::init_");
  }
  
}

// Remove references to the contained python objects
void PythonInstance::clean_()
{
  Py_DECREF(pLocals_);
  Py_DECREF(pCode_);
  
  if(Py_IsInitialized())
  {
    PyGC_Collect();
  }
}

// Call the python function object and return the result
PyObject* PythonInstance::call(PyObject *pArgs) const
{
  return PyObject_CallObject(pFunc_, pArgs);
}

