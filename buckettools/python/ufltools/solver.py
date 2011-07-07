from ufltools.base import *
import sys
import subprocess

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
    self.form_ranks = None
    self.system = None

  def ufl(self):
    """Write the system of forms to an array."""
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

    return ufl
    
  def namespace(self):
    return self.system.name+self.name

  def write_ufl(self, suffix=None):
    """Write the system of forms to a ufl file."""
    ufl = self.ufl()
    
    filename = self.namespace()+".ufl"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()
    
  def write_ufc(self, suffix=None):
    """Write the system of forms to a ufl file."""
    ufl = self.ufl()
    
    filename = self.namespace()+".ufl"
    if suffix: filename += "."+suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

    subprocess.call(["ffc", "-l", "dolfin", filename])
    
  def functionspace_cpp_no_if(self):
    return "      functionspace.reset( new "+self.namespace()+"::FunctionSpace(*mesh) );\n"

  def functionspace_cpp(self):
    cpp = [] 
    cpp.append("          case \""+self.name+"\":\n")
    cpp.append("            FunctionSpace_ptr functionspace(new "+self.system.name+self.name+"::FunctionSpace(*mesh));\n")
    cpp.append("            break;\n")
    return cpp

  def form_cpp(self):
    cpp = []
    for f in range(len(self.forms)):
      cpp.append("                  case \""+self.form_names[f]+"\":\n")
      if self.form_ranks[f]==0:
        cpp.append("                    Form_ptr form(new "+self.system.name+self.name+"::Form_"+`f`+"(*functionspace));\n")
      elif self.form_ranks[f]==1:
        cpp.append("                    Form_ptr form(new "+self.system.name+self.name+"::Form_"+`f`+"(*functionspace, *functionspace));\n")
      else:
        print "Unknwon form rank."
        sys.exit(1)
      cpp.append("                    break;\n")
    return cpp
      

