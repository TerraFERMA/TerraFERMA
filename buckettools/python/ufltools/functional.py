from ufltools.base import *

class Functional:
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl)."""

  def __init__(self):
    """Define the expected members of the functional class."""
    self.form = None
    self.function = None
    self.name = None
    self.symbol = None
  
  def ufl(self):
    """Write the functional to an array of ufl strings."""
    ufl = []
    if self.function.type=="Constant":
      ufl.append(declaration_comment("Coefficient", self.function.type, self.function.name))
      ufl.append(constant_ufl(self.function.symbol, self.function.system.cell))
    else:
      ufl += self.function.element_ufl()
      ufl.append("\n")
      ufl.append(declaration_comment("Coefficient", self.function.type, self.function.name))
      ufl.append(coefficient_ufl(self.function.symbol))
    ufl.append("\n")
    ufl.append(declaration_comment("Form", "form", self.name))
    ufl.append(self.form)
    ufl.append("\n")
    ufl.append(comment("Declare potentially non-default form names to be accessible"))
    ufl.append(forms_ufl([self.symbol]))
    ufl.append("\n")
    ufl.append(produced_comment())

    return ufl

  def write_ufl(self, suffix=None):
    """Write the functional to a ufl file."""
    ufl = self.ufl()

    filename   = self.function.system.name+self.function.name+self.name+".ufl"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

  def cpp(self):
    cpp = [] 
    cpp.append("              case \""+self.name+"\":\n")
    cpp.append("                Form_ptr functional(new "+self.function.system.name+self.function.name+self.name+"::Form_0(*functionspace));\n")
    cpp.append("                break;\n")
    return cpp

