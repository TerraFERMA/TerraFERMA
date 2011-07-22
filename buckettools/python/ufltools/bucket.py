from ufltools.base import *

class Bucket:
  """A class that stores all the information necessary to write the ufl for an options file (i.e. set of mixed function spaces)."""

  def __init__(self):
    """Define the expected members of the bucket class - only one really."""
    self.systems = None

  def write_ufc(self):
    """Write all ufl files described by the bucket."""
    for system in self.systems:
      for field in system.fields:
        for functional in field.functionals:
          functional.write_ufc()
      for coeff in system.coeffs:
        for functional in coeff.functionals:
          functional.write_ufc()
      for solver in system.solvers:
        solver.write_ufc()

  def write_ufl(self):
    """Write all ufl files described by the bucket."""
    for system in self.systems:
      for field in system.fields:
        for functional in field.functionals:
          functional.write_ufl()
      for coeff in system.coeffs:
        for functional in coeff.functionals:
          functional.write_ufl()
      for solver in system.solvers:
        solver.write_ufl()

  def write_cpp(self):
    """Write a cpp header file describing all the namespaces in the bucket."""
 
    cpp       = []
    cpp.append("#ifndef __SYSTEMS_WRAPPER_H\n")
    cpp.append("#define __SYSTEMS_WRAPPER_H\n")
    cpp.append("\n")
    cpp.append("#include \"SystemsWrapper.h\"\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    functionspace_cpp         = []
    functionspace_cpp.append("  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).\n")
    functionspace_cpp.append("  FunctionSpace_ptr ufc_fetch_functionspace(std::string systemname, Mesh_ptr mesh)\n")
    functionspace_cpp.append("  {\n")
    functionspace_cpp.append("    FunctionSpace_ptr functionspace;\n")

    solverfunctionspace_cpp         = []
    solverfunctionspace_cpp.append("  // A function to return a functionspace from a system given a mesh and a solvername.\n")
    solverfunctionspace_cpp.append("  FunctionSpace_ptr ufc_fetch_functionspace(std::string systemname, std::string solvername, Mesh_ptr mesh)\n")
    solverfunctionspace_cpp.append("  {\n")
    solverfunctionspace_cpp.append("    FunctionSpace_ptr functionspace;\n")

    solvercoefficientspace_cpp         = []
    solvercoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.\n")
    solvercoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace(std::string systemname, std::string solvername, std::string coefficientname, Mesh_ptr mesh)\n")
    solvercoefficientspace_cpp.append("  {\n")
    solvercoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functionalcoefficientspace_cpp         = []
    functionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.\n")
    functionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace(std::string systemname, std::string functionname, std::string functionalname, std::string coefficientname, Mesh_ptr mesh)\n")
    functionalcoefficientspace_cpp.append("  {\n")
    functionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functional_cpp            = []
    functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
    functional_cpp.append("  Form_ptr ufc_fetch_functional(std::string systemname, std::string functionname, std::string functionalname, Mesh_ptr mesh)\n")
    functional_cpp.append("  {\n")
    functional_cpp.append("    Form_ptr functional;\n")

    form_cpp = []
    form_cpp.append("  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.\n")
    form_cpp.append("  Form_ptr ufc_fetch_form(std::string systemname, std::string solvername, std::string solvertype, std::string formname, FunctionSpace_ptr functionspace)\n")
    form_cpp.append("  {\n")
    form_cpp.append("    Form_ptr form;\n")
 
    for s in range(len(self.systems)):
      include_cpp    += self.systems[s].include_cpp()
      functionspace_cpp += self.systems[s].functionspace_cpp(index=s)
      solverfunctionspace_cpp += self.systems[s].solverfunctionspace_cpp(index=s)
      solvercoefficientspace_cpp += self.systems[s].solvercoefficientspace_cpp(index=s)
      functionalcoefficientspace_cpp += self.systems[s].functionalcoefficientspace_cpp(index=s)
      functional_cpp += self.systems[s].functional_cpp(index=s)
      form_cpp += self.systems[s].form_cpp(index=s)

    functionspace_cpp.append("    else\n")
    functionspace_cpp.append("    {\n")
    functionspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_functionspace\");\n")
    functionspace_cpp.append("    }\n")
    functionspace_cpp.append("    return functionspace;\n")
    functionspace_cpp.append("  }\n")

    solverfunctionspace_cpp.append("    else\n")
    solverfunctionspace_cpp.append("    {\n")
    solverfunctionspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_functionspace\");\n")
    solverfunctionspace_cpp.append("    }\n")
    solverfunctionspace_cpp.append("    return functionspace;\n")
    solverfunctionspace_cpp.append("  }\n")

    solvercoefficientspace_cpp.append("    else\n")
    solvercoefficientspace_cpp.append("    {\n")
    solvercoefficientspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_coefficientspace\");\n")
    solvercoefficientspace_cpp.append("    }\n")
    solvercoefficientspace_cpp.append("    return coefficientspace;\n")
    solvercoefficientspace_cpp.append("  }\n")

    functionalcoefficientspace_cpp.append("    else\n")
    functionalcoefficientspace_cpp.append("    {\n")
    functionalcoefficientspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_coefficientspace\");\n")
    functionalcoefficientspace_cpp.append("    }\n")
    functionalcoefficientspace_cpp.append("    return coefficientspace;\n")
    functionalcoefficientspace_cpp.append("  }\n")

    functional_cpp.append("    else\n")
    functional_cpp.append("    {\n")
    functional_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_functional\");\n")
    functional_cpp.append("    }\n")
    functional_cpp.append("    return functional;\n")
    functional_cpp.append("  }\n")

    form_cpp.append("    else\n")
    form_cpp.append("    {\n")
    form_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_form\");\n")
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
    cpp += solvercoefficientspace_cpp
    cpp.append("\n")
    cpp += functionalcoefficientspace_cpp
    cpp.append("\n")
    cpp += form_cpp
    cpp.append("\n")
    cpp += functional_cpp
    cpp.append("}\n")
    cpp.append("\n")
    cpp.append("#endif\n")

    filename = "SystemsWrapper.cpp"
    filehandle = file(filename, 'w')
    filehandle.writelines(cpp)
    filehandle.close()

