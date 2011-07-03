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
        ufl.append(constant_ufl(coeff.symbol, coeff.system.cell))
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

