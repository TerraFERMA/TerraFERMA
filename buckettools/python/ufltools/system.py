from ufltools.base import *

class System:
  def __init__(self):
    self.cell = None
    self.mesh_name = None
    self.name = None
    self.symbol = None
    self.fields = []
    self.coeffs = []

  def functions_ufl(self):
    ufl = []
    for field in self.fields:
      ufl += field.element_ufl()
    ufl += self.element_ufl()
    ufl += self.test_ufl()
    ufl += self.trial_ufl()
    ufl += self.iterate_ufl()
    ufl += self.old_ufl()
    for coeff in self.coeffs:
      if coeff.type == "Constant":
        ufl.append(constant_ufl(coeff.symbol, coeff.system.cell))
      else:
        ufl += coeff.element_ufl()
        ufl.append(coefficient_ufl(coeff.symbol))
    ufl.append("\n")
    return ufl

  def element_ufl(self):
    ufl = []
    if len(self.fields)==1:
      ufl.append(equal_ufl(self.symbol, self.fields[0].symbol, suffix="_e")+"\n")
    else:
      ufl = self.symbol+"_e = "
      ufl += "MixedElement(["
      for s in range(len(self.fields)-1):
        ufl += self.fields[s].symbol+"_e, "
      ufl += self.fields[-1].symbol+"_e"
      ufl += "])\n"
    return ufl

  def test_ufl(self):
    ufl = []
    if len(self.fields)==1:
      ufl.append(testfunction_ufl(self.symbol))
      ufl.append(testfunction_ufl(self.fields[0].symbol))
    else:
      ufl.append(testfunction_ufl(self.symbol))
      ufl.append(self.split_ufl(suffix="_t"))
    return ufl

  def trial_ufl(self):
    ufl = []
    if len(self.fields)==1:
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(trialfunction_ufl(self.fields[0].symbol))
    else:
      ufl.append(trialfunction_ufl(self.symbol))
      ufl.append(self.split_ufl(suffix="_a"))
    return ufl

  def iterate_ufl(self):
    ufl = []
    if len(self.fields)==1:
      ufl.append(coefficient_ufl(self.symbol, suffix="_i"))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_i"))
    else:
      ufl.append(coefficient_ufl(self.symbol, suffix="_i"))
      ufl.append(self.split_ufl(suffix="_i"))
    return ufl

  def old_ufl(self):
    ufl = []
    if len(self.fields)==1:
      ufl.append(coefficient_ufl(self.symbol, suffix="_n"))
      ufl.append(coefficient_ufl(self.fields[0].symbol, suffix="_n"))
    else:
      ufl.append(coefficient_ufl(self.symbol, suffix="_n"))
      ufl.append(self.split_ufl(suffix="_n"))
    return ufl

  def split_ufl(self, suffix=""):
    ufl = "("
    for s in range(len(self.fields)-1):
      ufl += self.fields[s].symbol+suffix+", "
    ufl += self.fields[-1].symbol+suffix
    ufl += ") = split("+self.symbol+suffix+")\n"
    return ufl

