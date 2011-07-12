
#include "PythonInstance.h"
#include "Python.h"
#include <string>
#include <dolfin.h>

using namespace buckettools;

PythonInstance::PythonInstance(std::string function) : function_(function)
{
  init_();
}

PythonInstance::~PythonInstance()
{
  clean_();
}

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

void PythonInstance::clean_()
{
  Py_DECREF(pLocals_);
  Py_DECREF(pCode_);
  
  if(Py_IsInitialized())
  {
    PyGC_Collect();
  }
}

std::string PythonInstance::function() const
{
  return function_;
}

PyObject* PythonInstance::call(PyObject *pArgs) const
{
  return PyObject_CallObject(pFunc_, pArgs);
}

