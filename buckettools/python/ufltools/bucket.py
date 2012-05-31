from ufltools.base import *
import shutil
import hashlib
import os
import sys
import subprocess

class Bucket:
  """A class that stores all the information necessary to write the ufl for an options file (i.e. set of mixed function spaces)."""

  def __init__(self):
    """Define the expected members of the bucket class - only one really."""
    self.parameters = None
    self.systems = None

  def write_ufc(self):
    """Write all ufc files described by the bucket."""
    # Loop over the systems
    for system in self.systems:
      for field in system.fields:
        for functional in field.functionals:
          functional.write_ufc()
      for coeff in system.coeffs:
        if coeff.functional:
          coeff.functional.write_ufc()
        for functional in coeff.functionals:
          functional.write_ufc()
      for solver in system.solvers:
        solver.write_ufc()

  def write_cppexpressions(self):
    """Write all cpp expression header files described by the bucket."""
    for system in self.systems:
      for field in system.fields:
        for cppexpression in field.cpp:
          cppexpression.write_cppheader_md5()
      for coeff in system.coeffs:
        for cppexpression in coeff.cpp:
          cppexpression.write_cppheader_md5()

  def write_ufl(self):
    """Write all ufl files described by the bucket."""
    # Write simple ufc files for visualization functionspaces
    # Loop over the systems
    for system in self.systems:
      for field in system.fields:
        for functional in field.functionals:
          functional.write_ufl()
      for coeff in system.coeffs:
        if coeff.functional:
          coeff.functional.write_ufl()
        for functional in coeff.functionals:
          functional.write_ufl()
      for solver in system.solvers:
        solver.write_ufl()

  def write_systemfunctionals_cpp(self):
    """Write a cpp header file describing all the ufc namespaces in the bucket."""
    cpp = []

    cpp.append("\n")
    cpp.append("#include \"SystemFunctionalsWrapper.h\"\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    functionalcoefficientspace_cpp         = []
    functionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a functionname and a uflsymbol.\n")
    functionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    functionalcoefficientspace_cpp.append("  {\n")
    functionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    constantfunctionalcoefficientspace_cpp         = []
    constantfunctionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a coefficientname and a uflsymbol.\n")
    constantfunctionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &coefficientname, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    constantfunctionalcoefficientspace_cpp.append("  {\n")
    constantfunctionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functional_cpp            = []
    functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
    functional_cpp.append("  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, Mesh_ptr mesh)\n")
    functional_cpp.append("  {\n")
    functional_cpp.append("    Form_ptr functional;\n")

    constantfunctional_cpp            = []
    constantfunctional_cpp.append("  // A function to return a functional for a constant from a system-function set given a mesh.\n")
    constantfunctional_cpp.append("  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &coefficientname, Mesh_ptr mesh)\n")
    constantfunctional_cpp.append("  {\n")
    constantfunctional_cpp.append("    Form_ptr functional;\n")

    s = 0
    for system in self.systems:
      include_cpp    += system.include_systemfunctionals_cpp()
      functionalcoefficientspace_cpp += system.functionalcoefficientspace_cpp(index=s)
      constantfunctionalcoefficientspace_cpp += system.constantfunctionalcoefficientspace_cpp(index=s)
      functional_cpp += system.functional_cpp(index=s)
      constantfunctional_cpp += system.constantfunctional_cpp(index=s)
      s += 1

    functionalcoefficientspace_cpp.append("    else\n")
    functionalcoefficientspace_cpp.append("    {\n")
    functionalcoefficientspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_coefficientspace_from_functional\");\n")
    functionalcoefficientspace_cpp.append("    }\n")
    functionalcoefficientspace_cpp.append("    return coefficientspace;\n")
    functionalcoefficientspace_cpp.append("  }\n")

    constantfunctionalcoefficientspace_cpp.append("    else\n")
    constantfunctionalcoefficientspace_cpp.append("    {\n")
    constantfunctionalcoefficientspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_coefficientspace_from_functional\");\n")
    constantfunctionalcoefficientspace_cpp.append("    }\n")
    constantfunctionalcoefficientspace_cpp.append("    return coefficientspace;\n")
    constantfunctionalcoefficientspace_cpp.append("  }\n")

    functional_cpp.append("    else\n")
    functional_cpp.append("    {\n")
    functional_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_functional\");\n")
    functional_cpp.append("    }\n")
    functional_cpp.append("    return functional;\n")
    functional_cpp.append("  }\n")

    constantfunctional_cpp.append("    else\n")
    constantfunctional_cpp.append("    {\n")
    constantfunctional_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_functional\");\n")
    constantfunctional_cpp.append("    }\n")
    constantfunctional_cpp.append("    return functional;\n")
    constantfunctional_cpp.append("  }\n")

    cpp += include_cpp
    cpp.append("\n")
    cpp.append("namespace buckettools\n")
    cpp.append("{\n")
    cpp += functionalcoefficientspace_cpp
    cpp.append("\n")
    cpp += constantfunctionalcoefficientspace_cpp
    cpp.append("\n")
    cpp += functional_cpp
    cpp.append("\n")
    cpp += constantfunctional_cpp
    cpp.append("\n")
    cpp.append("}\n")
    cpp.append("\n")

    filename = "SystemFunctionalsWrapper.cpp"
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

  def write_systemsolvers_cpp(self):
    """Write a cpp header file describing all the ufc namespaces in the bucket."""
    cpp = []

    cpp.append("\n")
    cpp.append("#include \"SystemSolversWrapper.h\"\n")
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
    solvercoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace_from_solver(const std::string &systemname, const std::string &solvername, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    solvercoefficientspace_cpp.append("  {\n")
    solvercoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    form_cpp = []
    form_cpp.append("  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.\n")
    form_cpp.append("  Form_ptr ufc_fetch_form(const std::string &systemname, const std::string &solvername, const std::string &solvertype, const std::string &formname, const FunctionSpace_ptr functionspace)\n")
    form_cpp.append("  {\n")
    form_cpp.append("    Form_ptr form;\n")
 
    s = 0
    for system in self.systems:
      include_cpp    += system.include_systemsolvers_cpp()
      functionspace_cpp += system.functionspace_cpp(index=s)
      solverfunctionspace_cpp += system.solverfunctionspace_cpp(index=s)
      solvercoefficientspace_cpp += system.solvercoefficientspace_cpp(index=s)
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
    solvercoefficientspace_cpp.append("      dolfin::error(\"Unknown systemname in ufc_fetch_coefficientspace_from_solver\");\n")
    solvercoefficientspace_cpp.append("    }\n")
    solvercoefficientspace_cpp.append("    return coefficientspace;\n")
    solvercoefficientspace_cpp.append("  }\n")

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
    cpp += form_cpp
    cpp.append("\n")
    cpp.append("}\n")
    cpp.append("\n")

    filename = "SystemSolversWrapper.cpp"
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

  def write_systemexpressions_cpp(self):
    """Write a cpp header file describing all the cpp expression namespaces in the bucket."""
    cpp = []

    cpp.append("\n")
    cpp.append("#include \"SystemExpressionsWrapper.h\"\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    cppexpression_cpp = []
    cppexpression_cpp.append("  // A function to return an expression for a coefficient from a system given a systemname and a functionname (and its size, shape and private members bucket, system and time.\n")
    cppexpression_cpp.append("  Expression_ptr cpp_fetch_expression(const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname, const uint &size, const std::vector<int> &shape, const Bucket *bucket, const SystemBucket *system, const double_ptr time)\n")
    cppexpression_cpp.append("  {\n")
    cppexpression_cpp.append("    Expression_ptr expression;\n")
 
    cppexpression_init = []
    cppexpression_init.append("  // A function to initialize an expression for a cpp expression given a systemname and a functionname (and a boost shared pointer to the expression to initialize.\n")
    cppexpression_init.append("  void cpp_init_expression(Expression_ptr expression, const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname)\n")
    cppexpression_init.append("  {\n")
 
    s = 0
    for system in self.systems:
      include_cpp    += system.include_systemexpressions_cpp()
      cppexpression_cpp += system.cppexpression_cpp(index=s)
      cppexpression_init += system.cppexpression_init(index=s)
      s += 1

    cppexpression_cpp.append("    else\n")
    cppexpression_cpp.append("    {\n")
    cppexpression_cpp.append("      dolfin::error(\"Unknown systemname in cpp_fetch_expression\");\n")
    cppexpression_cpp.append("    }\n")
    cppexpression_cpp.append("    return expression;\n")
    cppexpression_cpp.append("  }\n")

    cppexpression_init.append("    else\n")
    cppexpression_init.append("    {\n")
    cppexpression_init.append("      dolfin::error(\"Unknown systemname in cpp_init_expression\");\n")
    cppexpression_init.append("    }\n")
    cppexpression_init.append("  }\n")

    cpp += include_cpp
    cpp.append("\n")
    cpp.append("namespace buckettools\n")
    cpp.append("{\n")
    cpp += cppexpression_cpp
    cpp.append("\n")
    cpp += cppexpression_init
    cpp.append("\n")
    cpp.append("}\n")
    cpp.append("\n")

    filename = "SystemExpressionsWrapper.cpp"
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


