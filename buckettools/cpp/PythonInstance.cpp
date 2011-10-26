
#include "PythonInstance.h"
#include "Python.h"
#include <string>
#include <dolfin.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PythonInstance::PythonInstance(const std::string &function) : 
                                                  function_(function)
{
  init_();                                                           // initialize
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PythonInstance::~PythonInstance()
{
  clean_();                                                          // finalize
}

//*******************************************************************|************************************************************//
// given a python arguments object call the val function and return the result
//*******************************************************************|************************************************************//
PyObject* PythonInstance::call(PyObject *pArgs) const
{
  return PyObject_CallObject(pFunc_, pArgs);
}

//*******************************************************************|************************************************************//
// initialize the python objects contained in the python instance
//*******************************************************************|************************************************************//
void PythonInstance::init_()
{
  if(!Py_IsInitialized())
  {
    Py_Initialize();                                                 // if python isn't initialized, do it now
  }
  
  pMain_ = PyImport_AddModule("__main__");                           // set up a main module
  pGlobals_ = PyModule_GetDict(pMain_);                              // set up the globals dictionary
  pLocals_ = PyDict_New();                                           // set up the locals dictionary
  
  pCode_ = PyRun_String((char*)function_.c_str(), Py_file_input,     // run the function string 
                                              pGlobals_, pLocals_);
  
  pFunc_ = PyDict_GetItemString(pLocals_, "val");                    // get the val function from the function string

  std::stringstream pythonbuffer;                                    // set up a simple python command to check how many 
  pythonbuffer << "import inspect" << std::endl                      // arguments the python function val has
               << "_nargs = len(inspect.getargspec(val).args)" 
               << std::endl;
  PyObject* tmppCode = PyRun_String(pythonbuffer.str().c_str(),      // run the python commands
                                Py_file_input, pGlobals_, pLocals_); 
  PyObject* pnArgs = PyDict_GetItemString(pLocals_, "_nargs");       // retrieve the result, _nargs
  nargs_ = PyInt_AsLong(pnArgs);                                     // recast it as an integer
  
  if (PyErr_Occurred()){                                             // check for errors in getting the function
    PyErr_Print();
    dolfin::error("In PythonInstance::init_");
  }
  
}

//*******************************************************************|************************************************************//
// finalize the python instance
//*******************************************************************|************************************************************//
void PythonInstance::clean_()
{
  Py_DECREF(pLocals_);                                               // decrease the reference count on the locals
  Py_DECREF(pCode_);                                                 // and the code
  
  if(Py_IsInitialized())                                             // if python is initialized then collect (safe for multiple
  {                                                                  // python instances?)
    PyGC_Collect();
  }
}

