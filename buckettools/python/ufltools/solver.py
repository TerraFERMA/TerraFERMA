from ufltools.base import *

class Solver:
  """A class that stores all the information necessary to write the ufl for a system of forms (i.e. linear or bilinear) associated with a solver."""

  def __init__(self):
    """Define the expected members of the solver class."""
    self.name = None
    self.type = None
    self.preamble = None
    self.forms = None
    self.form_symbols = None
    self.form_names = None
    self.system = None

  def write_ufl(self):
    """Write the system of forms to a ufl file."""
    ufl = []
    ufl += self.system.functions_ufl()
    if self.preamble:
      ufl.append(comment("Form preamble"))
      ufl.append(self.preamble+"\n")

    assert(len(self.forms)==len(self.form_names))
    for i in range(len(self.forms)):
      ufl.append(declaration_comment("Form", "form", self.form_names[i]))
      ufl.append(self.forms[i]+"\n")

    ufl.append("\n")
    assert(len(self.forms)==len(self.form_symbols))
    ufl.append(comment("Declare potentially non-default form names to be accessible"))
    ufl.append(forms_ufl(self.form_symbols))
    ufl.append("\n")
    ufl.append(produced_comment())

    filename = self.system.name+self.name+".ufl"
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()
    
