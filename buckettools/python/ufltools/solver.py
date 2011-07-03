from ufltools.base import *

class Solver:
  def __init__(self):
    self.name = None
    self.type = None
    self.preamble = None
    self.forms = None
    self.form_symbols = None
    self.form_names = None
    self.system = None

  def write_ufl(self):
    ufl = []
    ufl += self.system.functions_ufl()
    if self.preamble:
      ufl.append(generic_comment("Form preamble"))
      ufl.append(self.preamble+"\n")

    assert(len(self.forms)==len(self.form_names))
    for i in range(len(self.forms)):
      ufl.append(declaration_comment("Form", self.form_names[i]))
      ufl.append(self.forms[i]+"\n")

    ufl.append("\n")
    assert(len(self.forms)==len(self.form_symbols))
    ufl.append(forms_ufl(self.form_symbols))
    ufl.append("\n")
    ufl.append(produced_comment())

    filename = self.system.name+self.name+".ufl"
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()
    
