from ufltools.base import *
import shutil
import hashlib

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
    functionspace_cpp.append("  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname, Mesh_ptr mesh)\n")
    functionspace_cpp.append("  {\n")
    functionspace_cpp.append("    FunctionSpace_ptr functionspace;\n")

    solverfunctionspace_cpp         = []
    solverfunctionspace_cpp.append("  // A function to return a functionspace from a system given a mesh and a solvername.\n")
    solverfunctionspace_cpp.append("  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname, const std::string &solvername, Mesh_ptr mesh)\n")
    solverfunctionspace_cpp.append("  {\n")
    solverfunctionspace_cpp.append("    FunctionSpace_ptr functionspace;\n")

    solvercoefficientspace_cpp         = []
    solvercoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a uflsymbol.\n")
    solvercoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace(const std::string &systemname, const std::string &solvername, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    solvercoefficientspace_cpp.append("  {\n")
    solvercoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functionalcoefficientspace_cpp         = []
    functionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a uflsymbol.\n")
    functionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace(const std::string &systemname, const std::string &uflsymbol, const std::string &functionalname, const std::string &coefficientname, Mesh_ptr mesh)\n")
    functionalcoefficientspace_cpp.append("  {\n")
    functionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functional_cpp            = []
    functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
    functional_cpp.append("  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, Mesh_ptr mesh)\n")
    functional_cpp.append("  {\n")
    functional_cpp.append("    Form_ptr functional;\n")

    form_cpp = []
    form_cpp.append("  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.\n")
    form_cpp.append("  Form_ptr ufc_fetch_form(const std::string &systemname, const std::string &solvername, const std::string &solvertype, const std::string &formname, const FunctionSpace_ptr functionspace)\n")
    form_cpp.append("  {\n")
    form_cpp.append("    Form_ptr form;\n")
 
    s = 0
    for system in self.systems:
      include_cpp    += system.include_cpp()
      functionspace_cpp += system.functionspace_cpp(index=s)
      solverfunctionspace_cpp += system.solverfunctionspace_cpp(index=s)
      solvercoefficientspace_cpp += system.solvercoefficientspace_cpp(index=s)
      functionalcoefficientspace_cpp += system.functionalcoefficientspace_cpp(index=s)
      functional_cpp += system.functional_cpp(index=s)
      form_cpp += system.form_cpp(index=s)
      s += 1

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
    filehandle = file(filename+".temp", 'w')
    filehandle.writelines(cpp)
    filehandle.close()

    try:
      checksum = hashlib.md5(open(filename).read()).hexdigest()
    except:
      checksum = None

    if checksum != hashlib.md5(open(filename+".temp").read()).hexdigest():
      # file has changed
      shutil.copy(filename+".temp", filename)
