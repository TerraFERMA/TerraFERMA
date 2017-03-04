# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

from buckettools.base import *
import shutil
import hashlib
import os
import sys
import subprocess

class Bucket:
  """A class that stores all the information necessary to write the ufl for an options file (i.e. set of mixed function spaces)."""

  def __init__(self):
    """Define the expected members of the bucket class - only one really."""
    self.meshes = None
    self.viselementfamily = None
    self.viselementdegree = None
    self.parameters = None
    self.systems = None

  def write_ufc(self):
    """Write all ufc files described by the bucket."""
    # Write simple ufc files for visualization functionspaces
    self.write_visualization_ufcs()
    # Loop over the systems
    for system in self.systems:
      for coeff in system.coeffs:
        if coeff.functional:
          coeff.functional.write_ufc()
      for solver in system.solvers:
        solver.write_ufc()
      for functional in system.functionals:
        functional.write_ufc()

  def list_namespaces(self):
    """Return a list of the namespaces."""
    namespaces = []
    for meshname in self.meshes.iterkeys():
      namespaces.append(self.visualization_namespace(meshname))
    for system in self.systems:
      for coeff in system.coeffs:
        if coeff.functional:
          namespaces.append(coeff.functional.namespace())
      for solver in system.solvers:
        namespaces.append(solver.namespace())
      for functional in system.functionals:
        namespaces.append(functional.namespace())
    return namespaces

  def list_globaluflsymbols(self):
    """Return a list of global ufl_symbols."""
    uflsymbols = []
    if len(self.systems) > 0:
      for coeff in self.systems[0].special_coeffs:
        uflsymbols.append(coeff.symbol)
    for system in self.systems:
      uflsymbols.append(system.symbol)
      for field in system.fields:
        uflsymbols.append(field.symbol)
      for coeff in system.coeffs:
        uflsymbols.append(coeff.symbol)
    return uflsymbols

  def preprocess_checks(self):
    """Run preprocessing checks."""
    stat = self.repeated_uflsymbol_check()
    return stat

  def repeated_uflsymbol_check(self):
    """Check for repeated ufl symbols."""
    stat = 0
    uflsymbols = self.list_globaluflsymbols()
    repeated_uflsymbols = set([s for s in uflsymbols if uflsymbols.count(s) > 1])
    if len(repeated_uflsymbols) > 0: stat = 1
    for s in repeated_uflsymbols: print "ERROR global ufl_symbol %s repeated! Change one of its instances."%(s)
    repeated_auto_uflsymbols = set([(s, s+a) for s in uflsymbols for a in uflsymbol_suffixes() if uflsymbols.count(s+a) >= 1 and a != ''])
    if len(repeated_auto_uflsymbols) > 0: stat = 1
    for s in repeated_auto_uflsymbols: print "ERROR ufl_symbol generated from global ufl_symbol %s conflicts with global ufl_symbol %s! Change global ufl_symbol %s to avoid reserved endings."%(s[0], s[1], s[1])
    for system in self.systems:
      for coeff in system.coeffs:
        if coeff.functional:
          stat = max(stat, coeff.functional.repeated_uflsymbol_check())
      for solver in system.solvers:
        stat = max(stat, solver.repeated_uflsymbol_check())
      for functional in system.functionals:
        stat = max(stat, functional.repeated_uflsymbol_check())
    return stat

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
    self.write_visualization_ufls()
    # Loop over the systems
    for system in self.systems:
      for coeff in system.coeffs:
        if coeff.functional:
          coeff.functional.write_ufl()
      for solver in system.solvers:
        solver.write_ufl()
      for functional in system.functionals:
        functional.write_ufl()

  def visualization_ufl(self, meshname, meshcell):
    """Return the ufl describing the visualization functionspace."""
    ufl = []
    ufl.append(declaration_comment("Element", "Function", "VisualizationOnMesh"+meshname))
    ufl.append("vis_e = FiniteElement(\""+self.viselementfamily+"\", " \
               +meshcell+", " \
               +`self.viselementdegree`+")\n")
    ufl.append("\n")
    ufl.append(declaration_comment("Test space", "Function", "VisualizationOnMesh"+meshname))
    ufl.append(testfunction_ufl("vis"))
    ufl.append(declaration_comment("Trial space", "Function", "VisualizationOnMesh"+meshname))
    ufl.append(trialfunction_ufl("vis"))
    ufl.append(declaration_comment("Form", "form", "Bilinear"))
    ufl.append("a = vis_t*vis_a*dx\n")
    ufl.append("forms = [a]\n")
    ufl.append("\n")
    ufl.append(produced_comment())
    
    return ufl
      
  def visualization_namespace(self, meshname):
    return "_VisualizationOnMesh"+meshname
    
  def write_visualization_ufl(self, meshname, meshcell, suffix=None):
    """Write the mesh visualization functionspace to a ufl file."""
    ufl = self.visualization_ufl(meshname, meshcell)
    
    filename = self.visualization_namespace(meshname)+".ufl"
    if suffix: filename += suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()
    
  def write_visualization_ufls(self, suffix=None):
    """Write the mesh visualization functionspaces to a ufl file."""
    # loop over all the meshes we recorded information about in case they have different cells
    for meshname, meshcell in self.meshes.iteritems():
      self.write_visualization_ufl(meshname, meshcell, suffix=suffix)
    
  def write_visualization_ufcs(self):
    """Write the mesh visualization functionspaces to a ufl file and transform it into ufc."""
    # loop over all the meshes we recorded information about in case they have different cells
    for meshname, meshcell in self.meshes.iteritems():
      self.write_visualization_ufl(meshname, meshcell, suffix=".temp")
      ffc(self.visualization_namespace(meshname), 'quadrature', 'default', None)

  def write_systemfunctionals_cpp(self):
    """Write a cpp header file describing all the ufc namespaces in the bucket."""
    cpp = []

    cpp.append("\n")
    cpp.append("#include \"SystemFunctionalsWrapper.h\"\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include \"Logger.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    functionalcoefficientspace_cpp         = []
    functionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a functionname and a uflsymbol.\n")
    functionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &functionalname, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    functionalcoefficientspace_cpp.append("  {\n")
    functionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    constantfunctionalcoefficientspace_cpp         = []
    constantfunctionalcoefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a coefficientname and a uflsymbol.\n")
    constantfunctionalcoefficientspace_cpp.append("  FunctionSpace_ptr ufc_fetch_coefficientspace_from_constant_functional(const std::string &systemname, const std::string &coefficientname, const std::string &uflsymbol, Mesh_ptr mesh)\n")
    constantfunctionalcoefficientspace_cpp.append("  {\n")
    constantfunctionalcoefficientspace_cpp.append("    FunctionSpace_ptr coefficientspace;\n")

    functional_cpp            = []
    functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
    functional_cpp.append("  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &functionalname, Mesh_ptr mesh)\n")
    functional_cpp.append("  {\n")
    functional_cpp.append("    Form_ptr functional;\n")

    constantfunctional_cpp            = []
    constantfunctional_cpp.append("  // A function to return a functional for a constant from a system-function set given a mesh.\n")
    constantfunctional_cpp.append("  Form_ptr ufc_fetch_constant_functional(const std::string &systemname, const std::string &coefficientname, Mesh_ptr mesh)\n")
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
    functionalcoefficientspace_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_coefficientspace_from_functional\", \"System name: %s\", systemname.c_str());\n")
    functionalcoefficientspace_cpp.append("    }\n")
    functionalcoefficientspace_cpp.append("    return coefficientspace;\n")
    functionalcoefficientspace_cpp.append("  }\n")

    constantfunctionalcoefficientspace_cpp.append("    else\n")
    constantfunctionalcoefficientspace_cpp.append("    {\n")
    constantfunctionalcoefficientspace_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_coefficientspace_from_constant_functional\", \"System name: %s\", systemname.c_str());\n")
    constantfunctionalcoefficientspace_cpp.append("    }\n")
    constantfunctionalcoefficientspace_cpp.append("    return coefficientspace;\n")
    constantfunctionalcoefficientspace_cpp.append("  }\n")

    functional_cpp.append("    else\n")
    functional_cpp.append("    {\n")
    functional_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_functional\", \"System name: %s\", systemname.c_str());\n")
    functional_cpp.append("    }\n")
    functional_cpp.append("    return functional;\n")
    functional_cpp.append("  }\n")

    constantfunctional_cpp.append("    else\n")
    constantfunctional_cpp.append("    {\n")
    constantfunctional_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_constant_functional\", \"System name: %s\", systemname.c_str());\n")
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
    cpp.append("#include \"Logger.h\"\n")
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
    functionspace_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_functionspace\", \"System name: %s\", systemname.c_str());\n")
    functionspace_cpp.append("    }\n")
    functionspace_cpp.append("    return functionspace;\n")
    functionspace_cpp.append("  }\n")

    solverfunctionspace_cpp.append("    else\n")
    solverfunctionspace_cpp.append("    {\n")
    solverfunctionspace_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_functionspace\", \"System name: %s\", systemname.c_str());\n")
    solverfunctionspace_cpp.append("    }\n")
    solverfunctionspace_cpp.append("    return functionspace;\n")
    solverfunctionspace_cpp.append("  }\n")

    solvercoefficientspace_cpp.append("    else\n")
    solvercoefficientspace_cpp.append("    {\n")
    solvercoefficientspace_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_coefficientspace_from_solver\", \"System name: %s\", systemname.c_str());\n")
    solvercoefficientspace_cpp.append("    }\n")
    solvercoefficientspace_cpp.append("    return coefficientspace;\n")
    solvercoefficientspace_cpp.append("  }\n")

    form_cpp.append("    else\n")
    form_cpp.append("    {\n")
    form_cpp.append("      tf_err(\"Unknown systemname in ufc_fetch_form\", \"System name: %s\", systemname.c_str());\n")
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

  def write_visualization_cpp(self):
    """Write a cpp header file describing all the ufc namespaces used for visualization in the bucket."""
    cpp = []

    cpp.append("\n")
    cpp.append("#include \"VisualizationWrapper.h\"\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    functionspace_cpp         = []
    functionspace_cpp.append("  // A function to return a functionspace for visualization given a mesh and a mesh name.\n")
    functionspace_cpp.append("  FunctionSpace_ptr ufc_fetch_visualization_functionspace(const std::string &meshname, Mesh_ptr mesh)\n")
    functionspace_cpp.append("  {\n")
    functionspace_cpp.append("    FunctionSpace_ptr functionspace;\n")

    s = 0
    for meshname in self.meshes.iterkeys():
      include_cpp.append("#include \""+self.visualization_namespace(meshname)+".h\"\n")
      
      if s == 0:
        functionspace_cpp.append("    if (meshname ==  \""+meshname+"\")\n")
      else:
        functionspace_cpp.append("    else if (meshname ==  \""+meshname+"\")\n")
      functionspace_cpp.append("    {\n")
      functionspace_cpp.append("      functionspace.reset( new "+self.visualization_namespace(meshname)+"::FunctionSpace(mesh) );\n")
      functionspace_cpp.append("    }\n")
      
      s += 1

    functionspace_cpp.append("    else\n")
    functionspace_cpp.append("    {\n")
    functionspace_cpp.append("      dolfin::error(\"Unknown meshname in ufc_fetch_visualization_functionspace\");\n")
    functionspace_cpp.append("    }\n")
    functionspace_cpp.append("    return functionspace;\n")
    functionspace_cpp.append("  }\n")

    cpp += include_cpp
    cpp.append("\n")
    cpp.append("namespace buckettools\n")
    cpp.append("{\n")
    cpp += functionspace_cpp
    cpp.append("\n")
    cpp.append("}\n")
    cpp.append("\n")

    filename = "VisualizationWrapper.cpp"
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
    cpp.append("#include \"Logger.h\"\n")
    cpp.append("#include <dolfin.h>\n")

    include_cpp = []

    cppexpression_cpp = []
    cppexpression_cpp.append("  // A function to return an expression for a coefficient from a system given a systemname and a functionname (and its size, shape and private members bucket, system and time.\n")
    cppexpression_cpp.append("  Expression_ptr cpp_fetch_expression(const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname, const std::size_t &size, const std::vector<std::size_t> &shape, const Bucket *bucket, const SystemBucket *system, const double_ptr time)\n")
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
    cppexpression_cpp.append("      tf_err(\"Unknown systemname in cpp_fetch_expression\", \"System name: %s\", systemname.c_str());\n")
    cppexpression_cpp.append("    }\n")
    cppexpression_cpp.append("    return expression;\n")
    cppexpression_cpp.append("  }\n")

    cppexpression_init.append("    else\n")
    cppexpression_init.append("    {\n")
    cppexpression_init.append("      tf_err(\"Unknown systemname in cpp_init_expression\", \"System name: %s\", systemname.c_str());\n")
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


