from ufltools.base import *

class Functional:
  def __init__(self):
    self.form = None
    self.function = None
    self.name = None
    self.symbol = None
  
  def write_ufl(self):
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

    filename   = self.function.system.name+self.function.name+self.name+".ufl"
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

