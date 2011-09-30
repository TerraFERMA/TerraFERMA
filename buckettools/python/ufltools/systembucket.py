from ufltools.base import *

class SystemBucket:
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space)."""

  def __init__(self):
    """Define the expected members of the system class."""
    self.cell = None
    self.mesh_name = None
    self.name = None
    self.symbol = None
    self.fields = None
    self.coeffs = None
    self.special_coeffs = None  # these are constants declared outside systems but still need system information
    self.solvers = None
    self.bucket = None

  def functions_ufl(self):
    """Write an array of ufl strings describing all the functions (fields and coefficients) within a system."""
    ufl = []
    ufl.append(declaration_comment("Function elements", "System", self.name))
    for field in self.fields:
      ufl += field.element_ufl()
    ufl += self.element_ufl()
    ufl.append("\n")
    ufl += self.test_ufl()
    ufl.append("\n")
    ufl += self.trial_ufl()
    ufl.append("\n")
    ufl += self.function_ufl()
    ufl.append("\n")
    ufl += self.iterate_ufl()
    ufl.append("\n")
    ufl += self.old_ufl()
    ufl.append("\n")
    for coeff in self.coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append("\n")
    # special_coeffs are only added for this system *not the other systems*
    ufl.append(comment("Declaring special coefficients, such as the timestep."))
    ufl.append("\n")
    for coeff in self.special_coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        if coeff.type == "Function":
          print "coefficient functionspaces not output for special coefficient functions"
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append("\n")
    ufl.append(comment("Finished declaring functions for this system, start on other systems."))
    ufl.append("\n")
    for system in self.bucket.systems:
      if system.name == self.name: continue
      ufl.append(declaration_comment("Function elements", "System", system.name))
      for field in system.fields:
        ufl += field.element_ufl()
      ufl += system.element_ufl()
      ufl.append("\n")
      ufl += system.function_ufl()
      ufl.append("\n")
      ufl += system.iterate_ufl()
      ufl.append("\n")
      ufl += system.old_ufl()
      ufl.append("\n")
      for coeff in system.coeffs:
        if coeff.type == "Constant":
          ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
          for suffix in uflsymbol_suffixes():
            ufl.append(coeff.constant_ufl(suffix=suffix))
        else:
          ufl += coeff.element_ufl()
          ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
          for suffix in uflsymbol_suffixes():
            ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append("\n")
    ufl.append(comment("Finished declaring functions for all other systems, start on forms."))
    ufl.append("\n")
    return ufl

  def element_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) element of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(comment("System element is not mixed"))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_e")+"\n")
    else:
      ufl.append(declaration_comment("Mixed element", "System", self.name))
      ufl_line = self.symbol+"_e = "
      ufl_line += "MixedElement(["
      for s in range(len(self.fields)-1):
        ufl_line += self.fields[s].symbol+"_e, "
      ufl_line += self.fields[-1].symbol+"_e"
      ufl_line += "])\n"
      ufl.append(ufl_line)
    return ufl

  def test_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) test space of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Test space", self.fields[0].type, self.fields[0].name))
      ufl.append(testfunction_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Test space", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_t")+"\n")
    else:
      ufl.append(declaration_comment("Test space", "System", self.name))
      ufl.append(testfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Test spaces", "functions", ""))
      ufl += self.split_ufl(suffix="_t")
    return ufl

  def trial_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) trial space of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Trial space", self.fields[0].type, self.fields[0].name))
      ufl.append(trialfunction_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_a")+"\n")
    else:
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Trial spaces", "functions", ""))
      ufl += self.split_ufl(suffix="_a")
    return ufl

  def function_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) field values of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol))
      ufl.append(declaration_comment("Value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol)+"\n")
    else:
      ufl.append(declaration_comment("Value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol))
      ufl.append(declaration_comment("Values", "functions", ""))
      ufl += self.split_ufl()
    return ufl

  def iterate_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) iterated field values of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Last iteration value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_i"))
      ufl.append(declaration_comment("Last iteration value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_i")+"\n")
    else:
      ufl.append(declaration_comment("Last iteration value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_i"))
      ufl.append(declaration_comment("Last iteration values", "functions", ""))
      ufl += self.split_ufl(suffix="_i")
    return ufl

  def old_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) old field values of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Previous time-level value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_n"))
      ufl.append(declaration_comment("Previous time-level value", "System", self.name))
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_n")+"\n")
    else:
      ufl.append(declaration_comment("Previous time-level value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_n"))
      ufl.append(declaration_comment("Previous time-level values", "functions", ""))
      ufl += self.split_ufl(suffix="_n")
    return ufl

  def split_ufl(self, suffix=""):
    """Write an array of ufl strings splitting a mixed function or element into components for its constitutive fields."""
    ufl = []
    for s in range(len(self.fields)):
      ufl.append(comment(" - "+self.fields[s].name))
    ufl_line = "("
    for s in range(len(self.fields)-1):
      ufl_line += self.fields[s].symbol+suffix+", "
    ufl_line += self.fields[-1].symbol+suffix
    ufl_line += ") = split("+self.symbol+suffix+")\n"
    ufl.append(ufl_line)
    return ufl

  def include_cpp(self):
    """Write an array of cpp strings including the namespaces of the system ufls."""
    cpp = []
    for field in self.fields:
      for functional in field.functionals:
        cpp.append("#include \""+functional.namespace()+".h\"\n")
    for coeff in self.coeffs:
      if coeff.functional:
        cpp.append("#include \""+coeff.functional.namespace()+".h\"\n")
      for functional in coeff.functionals:
        cpp.append("#include \""+functional.namespace()+".h\"\n")
    for solver in self.solvers:
      cpp.append("#include \""+solver.namespace()+".h\"\n")
    return cpp

  def functionspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functionspaces in the ufls (assuming a single solver)."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    cpp.append("      // All solvers within a system should return the same functionspace so just take the first one\n")
    cpp.append(self.solvers[0].functionspace_cpp_no_if())
    cpp.append("    }\n")
    return cpp

  def solverfunctionspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functionspaces in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    cpp.append("      // All solvers within a system should return the same functionspace\n")
    for s in range(len(self.solvers)):
      cpp += self.solvers[s].functionspace_cpp(index=s)
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        dolfin::error(\"Unknown solvername in ufc_fetch_functionspace\");\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    return cpp

  def functionalcoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    for c in range(len(self.fields)):
      if c == 0:
        cpp.append("      if (functionname ==  \""+self.fields[c].name+"\")\n")
      else:
        cpp.append("      else if (functionname ==  \""+self.fields[c].name+"\")\n")
      cpp.append("      {\n")
      if len(self.fields[c].functionals)==0:
        cpp.append("        dolfin::error(\"Unknown functionalname in ufc_fetch_coefficientspace_from_functional\");\n")
        cpp.append("      }\n")
      else:
        for f in range(len(self.fields[c].functionals)):
          if f == 0:
            cpp.append("        if (functionalname ==  \""+self.fields[c].functionals[f].name+"\")\n")
          else:
            cpp.append("        else if (functionalname ==  \""+self.fields[c].functionals[f].name+"\")\n")
          cpp.append("        {\n")
          symbols_in_form = form_symbols(self.fields[c].functionals[f].form)
          symbols_found = 0
          for c2 in range(len(self.coeffs)):
            for suffix in uflsymbol_suffixes():
              if self.coeffs[c2].symbol+suffix in symbols_in_form: 
                cpp += self.fields[c].functionals[f].coefficientspace_cpp(self.coeffs[c2], index=symbols_found, suffix=suffix)
                symbols_found =+ 1
                break
          if symbols_found == 0:
            cpp.append("          dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("        }\n")
          else:
            cpp.append("          else\n")
            cpp.append("          {\n")
            cpp.append("            dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("          }\n")
            cpp.append("        }\n")
        cpp.append("        else\n")
        cpp.append("        {\n")
        cpp.append("          dolfin::error(\"Unknown functionalname in ufc_fetch_coefficientspace_from_functional\");\n")
        cpp.append("        }\n")
        cpp.append("      }\n")
    for c in range(len(self.coeffs)):
      cpp.append("      else if (functionname ==  \""+self.coeffs[c].name+"\")\n")
      cpp.append("      {\n")
      if len(self.coeffs[c].functionals)==0:
        cpp.append("        dolfin::error(\"Unknown functionalname in ufc_fetch_coefficientspace_from_functional\");\n")
        cpp.append("      }\n")
      else:
        for f in range(len(self.coeffs[c].functionals)):
          if f == 0:
            cpp.append("        if (functionalname ==  \""+self.coeffs[c].functionals[f].name+"\")\n")
          else:
            cpp.append("        else if (functionalname ==  \""+self.coeffs[c].functionals[f].name+"\")\n")
          cpp.append("        {\n")
          symbols_in_form = form_symbols(self.coeffs[c].functionals[f].form)
          symbols_found = 0
          for c2 in range(len(self.coeffs)):
            for suffix in uflsymbol_suffixes():
              if self.coeffs[c2].symbol+suffix in symbols_in_form: 
                cpp += self.coeffs[c].functionals[f].coefficientspace_cpp(self.coeffs[c2], index=symbols_found, suffix=suffix)
                symbols_found =+ 1
                break
          if symbols_found == 0:
            cpp.append("          dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("        }\n")
          else:
            cpp.append("          else\n")
            cpp.append("          {\n")
            cpp.append("            dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("          }\n")
            cpp.append("        }\n")
        cpp.append("        else\n")
        cpp.append("        {\n")
        cpp.append("          dolfin::error(\"Unknown functionalname in ufc_fetch_coefficientspace_from_functional\");\n")
        cpp.append("        }\n")
        cpp.append("      }\n")
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        dolfin::error(\"Unknown functionname in ufc_fetch_coefficientspace_from_functional\");\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    return cpp

  def constantfunctionalcoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the constant functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    if (len(self.coeffs)==0):
      cpp.append("      dolfin::error(\"Unknown coefficentname in ufc_fetch_coefficientspace_from_functional\");\n")
      cpp.append("    }\n")
    else:
      for c in range(len(self.coeffs)):
        if c == 0:
          cpp.append("      if (coefficientname ==  \""+self.coeffs[c].name+"\")\n")
        else:
          cpp.append("      else if (coefficientname ==  \""+self.coeffs[c].name+"\")\n")
        cpp.append("      {\n")
        if self.coeffs[c].functional:
          symbols_in_form = form_symbols(self.coeffs[c].functional.form)
          symbols_found = 0
          for c2 in range(len(self.coeffs)):
            for suffix in uflsymbol_suffixes():
              if self.coeffs[c2].symbol+suffix in symbols_in_form: 
                cpp += self.coeffs[c].functional.coefficientspace_cpp(self.coeffs[c2], index=symbols_found, suffix=suffix)
                symbols_found =+ 1
                break
          if symbols_found == 0:
            cpp.append("        dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("      }\n")
          else:
            cpp.append("          else\n")
            cpp.append("          {\n")
            cpp.append("            dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional\");\n")
            cpp.append("          }\n")
            cpp.append("        }\n")
        else:
          cpp.append("        dolfin::error(\"Unknown functional in ufc_fetch_coefficientspace_from_functional\");\n")
          cpp.append("      }\n")
      cpp.append("      else\n")
      cpp.append("      {\n")
      cpp.append("        dolfin::error(\"Unknown coefficientname in ufc_fetch_coefficientspace_from_functional\");\n")
      cpp.append("      }\n")
      cpp.append("    }\n")
    return cpp

  def solvercoefficientspace_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    for s in range(len(self.solvers)):
      if s == 0:
        cpp.append("      if (solvername ==  \""+self.solvers[s].name+"\")\n")
      else:
        cpp.append("      else if (solvername ==  \""+self.solvers[s].name+"\")\n")
      cpp.append("      {\n")
      if len(self.coeffs)==0:
        cpp.append("        dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\");\n")
        cpp.append("      }\n")
      else:
        symbols_in_forms = forms_symbols(self.solvers[s].forms)
        if self.solvers[s].preamble: symbols_in_forms += form_symbols(self.solvers[s].preamble)
        symbols_found = 0
        for c in range(len(self.coeffs)):
          for suffix in uflsymbol_suffixes():
            if self.coeffs[c].symbol+suffix in symbols_in_forms: 
              cpp += self.coeffs[c].solvercoefficientspace_cpp(self.solvers[s].name, index=symbols_found, suffix=suffix)
              symbols_found =+ 1
              break
        if symbols_found == 0:
          cpp.append("        dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\");\n")
          cpp.append("      }\n")
        else:
          cpp.append("        else\n")
          cpp.append("        {\n")
          cpp.append("          dolfin::error(\"Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver\");\n")
          cpp.append("        }\n")
          cpp.append("      }\n")
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        dolfin::error(\"Unknown solvername in ufc_fetch_coefficientspace_from_solver\");\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    return cpp

  def form_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the forms in the ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    for s in range(len(self.solvers)):
      if s == 0:
        cpp.append("      if (solvername ==  \""+self.solvers[s].name+"\")\n")
      else:
        cpp.append("      else if (solvername ==  \""+self.solvers[s].name+"\")\n")
      cpp.append("      {\n")
      cpp.append("        if (solvertype == \""+self.solvers[s].type+"\")\n")
      cpp.append("        {\n")
      cpp += self.solvers[s].form_cpp()
      cpp.append("          else\n")
      cpp.append("          {\n")
      cpp.append("            dolfin::error(\"Unknown formname in ufc_fetch_form\");\n")
      cpp.append("          }\n")
      cpp.append("        }\n")
      cpp.append("        else\n")
      cpp.append("        {\n")
      cpp.append("          dolfin::error(\"Unknown solvertype in ufc_fetch_form\");\n")
      cpp.append("        }\n")
      cpp.append("      }\n")
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        dolfin::error(\"Unknown systemname in ufc_fetch_form\");\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    return cpp

  def functional_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    for f in range(len(self.fields)):
      if f == 0:
        cpp.append("      if (functionname ==  \""+self.fields[f].name+"\")\n")
      else:
        cpp.append("      else if (functionname ==  \""+self.fields[f].name+"\")\n")
      cpp.append("      {\n")
      for a in range(len(self.fields[f].functionals)):
        cpp += self.fields[f].functionals[a].cpp(index=a)
      if len(self.fields[f].functionals)==0:
        cpp.append("        dolfin::error(\"Unknown functionalname in ufc_fetch_functional\");\n")
        cpp.append("      }\n")
      else:
        cpp.append("        else\n")
        cpp.append("        {\n")
        cpp.append("          dolfin::error(\"Unknown functionalname in ufc_fetch_functional\");\n")
        cpp.append("        }\n")
        cpp.append("      }\n")
    for c in range(len(self.coeffs)):
      if c+len(self.fields) == 0:
        cpp.append("      if (functionname ==  \""+self.coeffs[c].name+"\")\n")
      else:
        cpp.append("      else if (functionname ==  \""+self.coeffs[c].name+"\")\n")
      cpp.append("      {\n")
      for a in range(len(self.coeffs[c].functionals)):
        cpp += self.coeffs[c].functionals[a].cpp(index=a)
      if len(self.coeffs[c].functionals)==0:
        cpp.append("        dolfin::error(\"Unknown functionalname in ufc_fetch_functional\");\n")
        cpp.append("      }\n")
      else:
        cpp.append("        else\n")
        cpp.append("        {\n")
        cpp.append("          dolfin::error(\"Unknown functionalname in ufc_fetch_functional\");\n")
        cpp.append("        }\n")
        cpp.append("      }\n")
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        dolfin::error(\"Unknown functionname in ufc_fetch_functional\");\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    return cpp

  def constantfunctional_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the constant coefficient functional ufls."""
    cpp = []  
    if index == 0:
      cpp.append("    if (systemname ==  \""+self.name+"\")\n")
    else:
      cpp.append("    else if (systemname ==  \""+self.name+"\")\n")
    cpp.append("    {\n")
    functionals_found = 0
    for c in range(len(self.coeffs)):
      if self.coeffs[c].functional:
        if functionals_found == 0:
          cpp.append("      if (coefficientname ==  \""+self.coeffs[c].name+"\")\n")
        else:
          cpp.append("      else if (coefficientname ==  \""+self.coeffs[c].name+"\")\n")
        functionals_found =+ 1
        cpp.append("      {\n")
        cpp.append("        functional.reset(new "+self.coeffs[c].functional.namespace()+"::Form_0(mesh));\n")
        cpp.append("      }\n")
    if functionals_found==0:
      cpp.append("      dolfin::error(\"Unknown coefficientname in ufc_fetch_functional\");\n")
      cpp.append("    }\n")
    else:
      cpp.append("      else\n")
      cpp.append("      {\n")
      cpp.append("        dolfin::error(\"Unknown coefficientname in ufc_fetch_functional\");\n")
      cpp.append("      }\n")
      cpp.append("    }\n")
    return cpp

