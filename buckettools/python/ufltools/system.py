from ufltools.base import *

class System:
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space)."""

  def __init__(self):
    """Define the expected members of the system class."""
    self.cell = None
    self.mesh_name = None
    self.name = None
    self.symbol = None
    self.fields = None
    self.coeffs = None
    self.solvers = None

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
    ufl += self.iterate_ufl()
    ufl.append("\n")
    ufl += self.old_ufl()
    ufl.append("\n")
    for coeff in self.coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        ufl.append(constant_ufl(coeff.symbol, self.cell))
      else:
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        ufl.append(coefficient_ufl(coeff.symbol))
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
      ufl.append(declaration_comment("Test space", "System", self.name))
      ufl.append(testfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Test space", self.fields[0].type, self.fields[0].name))
      ufl.append(testfunction_ufl(self.fields[0].symbol))
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
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Trial space", self.fields[0].type, self.fields[0].name))
      ufl.append(trialfunction_ufl(self.fields[0].symbol))
    else:
      ufl.append(declaration_comment("Trial space", "System", self.name))
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(declaration_comment("Trial spaces", "functions", ""))
      ufl += self.split_ufl(suffix="_a")
    return ufl

  def iterate_ufl(self):
    """Write an array of ufl strings describing the (potentially mixed) iterated field values of a system."""
    ufl = []
    if len(self.fields)==1:
      ufl.append(declaration_comment("Last iteration value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_i"))
      ufl.append(declaration_comment("Last iteration value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_i"))
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
      ufl.append(declaration_comment("Previous time-level value", "System", self.name))
      ufl.append(coefficient_ufl(self.symbol, suffix="_n"))
      ufl.append(declaration_comment("Previous time-level value", self.fields[0].type, self.fields[0].name))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_n"))
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
        cpp.append("#include \""+self.name+field.name+functional.name+".h\"\n")
    for solver in self.solvers:
      cpp.append("#include \""+self.name+solver.name+".h\"\n")
    return cpp

  def functionspace_cpp(self):
    """Write an array of cpp strings describing the namespace of the functionspaces in the ufls."""
    cpp = []  
    cpp.append("      case \""+self.name+"\":\n")
    cpp.append("        switch(solvername)\n")
    cpp.append("        {\n")
    for solver in self.solvers:
      cpp += solver.functionspace_cpp()
    cpp.append("          default:\n")
    cpp.append("            dolfin::error(\"Unknown solvername in fetch_functionspace\");\n")
    cpp.append("        }\n")
    cpp.append("        break;\n")
    return cpp

  def coefficientspace_cpp(self):
    """Write an array of cpp strings describing the namespace of the coefficientspaces in the ufls."""
    cpp = []  
    cpp.append("      case \""+self.name+"\":\n")
    cpp.append("        switch(solvername)\n")
    cpp.append("        {\n")
    for solver in self.solvers:
      cpp.append("        case \""+solver.name+"\":\n")
      cpp.append("          switch(functionname)\n")
      cpp.append("          {\n")
      for coeff in self.coeffs:
        if solver.preamble:
          coeff_present = any([coeff.symbol in form for form in solver.forms]+[coeff.symbol in solver.preamble])
        else:
          coeff_present = any([coeff.symbol in form for form in solver.forms])
        if coeff_present:
          cpp += coeff.coefficientspace_cpp(solver.name)
      cpp.append("            default:\n")
      cpp.append("              dolfin::error(\"Unknown functionname in fetch_coefficientspace\");\n")
      cpp.append("          }\n")
    cpp.append("          default:\n")
    cpp.append("            dolfin::error(\"Unknown solvername in fetch_coefficientspace\");\n")
    cpp.append("        }\n")
    cpp.append("        break;\n")
    return cpp

  def functional_cpp(self):
    """Write an array of cpp strings describing the namespace of the functional ufls."""
    cpp = []  
    cpp.append("      case \""+self.name+"\":\n")
    cpp.append("        switch(functionname)\n")
    cpp.append("        {\n")
    for field in self.fields:
      cpp.append("          case \""+field.name+"\":\n")
      cpp.append("            switch(functionalname)\n")
      cpp.append("            {\n")
      for functional in field.functionals:
        cpp += functional.cpp()
      cpp.append("              default:\n")
      cpp.append("                dolfin::error(\"Unknown functionalname in fetch_functional\");\n")
      cpp.append("            }\n")
      cpp.append("            break;\n")
    for coeff in self.coeffs:
      cpp.append("          case \""+coeff.name+"\":\n")
      cpp.append("            switch(functionalname)\n")
      cpp.append("            {\n")
      for functional in coeff.functionals:
        cpp += functional.cpp()
      cpp.append("              default:\n")
      cpp.append("                dolfin::error(\"Unknown functionalname in fetch_functional\");\n")
      cpp.append("            }\n")
      cpp.append("            break;\n")
    cpp.append("          default:\n")
    cpp.append("            dolfin::error(\"Unknown functionname in fetch_functional\");\n")
    cpp.append("        }\n")
    cpp.append("        break;\n")
    return cpp

