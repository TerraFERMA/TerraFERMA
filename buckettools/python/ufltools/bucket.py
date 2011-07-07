from ufltools.base import *

class Bucket:
  """A class that stores all the information necessary to write the ufl for an options file (i.e. set of mixed function spaces)."""

  def __init__(self):
    """Define the expected members of the bucket class - only one really."""
    self.systems = None

  def write_ufl(self, suffix=None):
    """Write all ufl files described by the bucket."""
    for system in self.systems:
      for field in system.fields:
        for functional in field.functionals:
          functional.write_ufl(suffix=suffix)
      for coeff in system.coeffs:
        for functional in coeff.functionals:
          functional.write_ufl(suffix=suffix)
      for solver in system.solvers:
        solver.write_ufl(suffix=suffix)

  def write_cpp(self, suffix=None):
    """Write a cpp header file describing all the namespaces in the bucket."""
 
    cpp       = []
    cpp.append("#ifndef __SYSTEMS_WRAPPER_H\n")
    cpp.append("#define __SYSTEMS_WRAPPER_H\n")
    cpp.append("\n")
    cpp.append("#include \"Bucket.h\"\n")
    cpp.append("#include <dolfin.h>\n")
    cpp.append("\n")

    include_cpp = []

    functionspace_cpp         = []
    functionspace_cpp.append("  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).\n")
    functionspace_cpp.append("  FunctionSpace_ptr fetch_functionspace(std::string systemname, Mesh_ptr mesh)\n")
    functionspace_cpp.append("  {\n")
    functionspace_cpp.append("    switch(systemname)\n")
    functionspace_cpp.append("    {\n")

    solverfunctionspace_cpp         = []
    solverfunctionspace_cpp.append("  // A function to return a functionspace from a system given a mesh and a solvername.\n")
    solverfunctionspace_cpp.append("  FunctionSpace_ptr fetch_functionspace(std::string systemname, std::string solvername, Mesh_ptr mesh)\n")
    solverfunctionspace_cpp.append("  {\n")
    solverfunctionspace_cpp.append("    switch(systemname)\n")
    solverfunctionspace_cpp.append("    {\n")

    coefficientspace_cpp         = []
    coefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.\n")
    coefficientspace_cpp.append("  FunctionSpace_ptr fetch_coefficientspace(std::string systemname, std::string solvername, std::string functionname, Mesh_ptr mesh)\n")
    coefficientspace_cpp.append("  {\n")
    coefficientspace_cpp.append("    switch(systemname)\n")
    coefficientspace_cpp.append("    {\n")

    functional_cpp            = []
    functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
    functional_cpp.append("  Form_ptr fetch_functional(std::string systemname, std::string functionname, std::string functionalname, Mesh_ptr mesh)\n")
    functional_cpp.append("  {\n")
    functional_cpp.append("    switch(systemname)\n")
    functional_cpp.append("    {\n")

    form_cpp = []
    form_cpp.append("  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.\n")
    form_cpp.append("  Form_ptr fetch_form(std::string systemname, std::string solvername, std::string solvertype, std::string formname, FunctionSpace_ptr functionspace)\n")
    form_cpp.append("  {\n")
    form_cpp.append("    switch(systemname)\n")
    form_cpp.append("    {\n")
 
    for system in self.systems:
      include_cpp    += system.include_cpp()
      functionspace_cpp += system.functionspace_cpp()
      solverfunctionspace_cpp += system.solverfunctionspace_cpp()
      coefficientspace_cpp += system.coefficientspace_cpp()
      functional_cpp += system.functional_cpp()
      form_cpp += system.form_cpp()

    functionspace_cpp.append("      default:\n")
    functionspace_cpp.append("        dolfin::error(\"Unknown systemname in fetch_functionspace\");\n")
    functionspace_cpp.append("    }\n")
    functionspace_cpp.append("    return functionspace;\n")
    functionspace_cpp.append("  }\n")

    solverfunctionspace_cpp.append("      default:\n")
    solverfunctionspace_cpp.append("        dolfin::error(\"Unknown systemname in fetch_functionspace\");\n")
    solverfunctionspace_cpp.append("    }\n")
    solverfunctionspace_cpp.append("    return functionspace;\n")
    solverfunctionspace_cpp.append("  }\n")

    coefficientspace_cpp.append("      default:\n")
    coefficientspace_cpp.append("        dolfin::error(\"Unknown systemname in fetch_coefficientspace\");\n")
    coefficientspace_cpp.append("    }\n")
    coefficientspace_cpp.append("    return coefficientspace;\n")
    coefficientspace_cpp.append("  }\n")

    functional_cpp.append("      default:\n")
    functional_cpp.append("        dolfin::error(\"Unknown systemname in fetch_functional\");\n")
    functional_cpp.append("    }\n")
    functional_cpp.append("    return functional;\n")
    functional_cpp.append("  }\n")

    form_cpp.append("      default:\n")
    form_cpp.append("        dolfin::error(\"Unknown systemname in fetch_form\");\n")
    form_cpp.append("    }\n")
    form_cpp.append("    return form;\n")
    form_cpp.append("  }\n")

    cpp += include_cpp
    cpp.append("\n")
    cpp.append("namespace buckettools\n")
    cpp.append("{\n")
    cpp += functionspace_cpp
    cpp.append("\n")
    cpp += solverfunctionspace_cpp
    cpp.append("\n")
    cpp += coefficientspace_cpp
    cpp.append("\n")
    cpp += form_cpp
    cpp.append("\n")
    cpp += functional_cpp
    cpp.append("}\n")
    cpp.append("\n")
    cpp.append("#endif\n")

    filename = "SystemsWrapper.h"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(cpp)
    filehandle.close()

