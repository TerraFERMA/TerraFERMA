from ufltools.base import *
import subprocess

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

    filename   = self.namespace()+".ufl"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

  def write_ufc(self, suffix=None):
    """Write the functional to a ufl file."""
    ufl = self.ufl()

    filename   = self.namespace()+".ufl"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

    subprocess.call(["ffc", "-l", "dolfin", filename])

  def namespace(self):
    return self.function.system.name+self.function.name+self.name

  def cpp(self, index=0):
    cpp = []
    if index == 0:
      cpp.append("        if (functionalname ==  \""+self.name+"\")\n")
    else:
      cpp.append("        else if (functionalname ==  \""+self.name+"\")\n")
    cpp.append("        {\n")
    cpp.append("          functional.reset(new "+self.namespace()+"::Form_0(*functionspace));\n")
    cpp.append("        }\n")
    return cpp

