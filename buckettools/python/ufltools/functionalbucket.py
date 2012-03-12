from ufltools.base import *
import subprocess
import hashlib
import shutil
import sys

class FunctionalBucket:
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl)."""

  def __init__(self):
    """Define the expected members of the functional class."""
    self.form = None
    self.function = None
    self.name = None
    self.symbol = None
    self.quadrature_degree = None
  
  def ufl(self):
    """Write the functional to an array of ufl strings."""
    ufl = []
    for field in self.function.system.fields:
      ufl += field.element_ufl()
    ufl += self.function.system.element_ufl()
    ufl.append("\n")
    ufl += self.function.system.function_ufl()
    ufl.append("\n")
    ufl += self.function.system.iterate_ufl()
    ufl.append("\n")
    ufl += self.function.system.old_ufl()
    ufl.append("\n")
    for coeff in self.function.system.coeffs:
      if coeff.type=="Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
      ufl.append("\n")
    ufl.append("\n")
    # special_coeffs are only added for this system *not the other systems*
    ufl.append(comment("Declaring special coefficients, such as the timestep."))
    ufl.append("\n")
    for coeff in self.function.system.special_coeffs:
      if coeff.type == "Constant":
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coeff.constant_ufl(suffix=suffix))
      else:
        if coeff.type == "Function":
          print "coefficient functionspaces not output for special coefficient functions"
          sys.exit(1)
        ufl += coeff.element_ufl()
        ufl.append(declaration_comment("Coefficient", coeff.type, coeff.name))
        for suffix in uflsymbol_suffixes():
          ufl.append(coefficient_ufl(coeff.symbol, suffix=suffix))
    ufl.append("\n")
    ufl.append(comment("Finished declaring functions for this system, start on other systems."))
    ufl.append("\n")
    for system in self.function.system.bucket.systems:
      if system.name == self.function.system.name: continue
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
    if self.function.system.bucket.parameters:
      ufl.append(comment("Global preamble"))
      ufl.append(self.function.system.bucket.parameters+"\n")
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
    if suffix: filename += suffix 
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()

  def write_ufc(self):
    """Write the functional to a ufl file."""
    self.write_ufl(suffix=".temp")

    filename = self.namespace()+".ufl"

    try:
      checksum = hashlib.md5(open(filename).read()).hexdigest()
    except:
      checksum = None

    try:
      headerfile = open(self.namespace()+".h")
    except IOError:
      headerfile = None

    if headerfile:
      rebuild = checksum != hashlib.md5(open(filename+".temp").read()).hexdigest()

      headertext = headerfile.read()
      qdi  = headertext.find("quadrature_degree")
      qdin = headertext.find("\n", qdi, -1)
      if (qdi != -1) and (qdin != -1):
        qd_old = headertext[qdi:qdin].split(" ")[-1]
        qd_new = "'"+`self.quadrature_degree`+"'" if self.quadrature_degree else "'auto'"
        rebuild = rebuild or qd_new != qd_old
      else: # if we don't know what the old quadrature degree was we always want to (re)build
        rebuild = True
    else:   # if we don't have a header file we always want to (re)build
      rebuild = True

    if rebuild:
      # files and/or quadrature_degree have changed
      shutil.copy(filename+".temp", filename)
      command = ["ffc", "-l", "dolfin", "-O", "-r", "quadrature"]
      if self.quadrature_degree:
        command += ["-f", "quadrature_degree="+`self.quadrature_degree`]
      command += [filename]
      try:
        subprocess.check_call(command)
      except:
        print "ERROR while calling ffc on file ", filename
        sys.exit(1)

  def namespace(self):
    return self.function.system.name+self.function.name+self.name

  def cpp(self, index=0):
    cpp = []
    if index == 0:
      cpp.append("        if (functionalname ==  \""+self.name+"\")\n")
    else:
      cpp.append("        else if (functionalname ==  \""+self.name+"\")\n")
    cpp.append("        {\n")
    cpp.append("          functional.reset(new "+self.namespace()+"::Form_0(mesh));\n")
    cpp.append("        }\n")
    return cpp

  def coefficientspace_cpp(self, coeff, index=0, suffix=""):
    cpp = [] 
    if index == 0:
      cpp.append("          if (uflsymbol ==  \""+coeff.symbol+"\")\n")
    else:
      cpp.append("          else if (uflsymbol ==  \""+coeff.symbol+"\")\n")
    cpp.append("          {\n")
    cpp.append("            coefficientspace.reset(new "+self.namespace()+"::CoefficientSpace_"+coeff.symbol+suffix+"(mesh));\n")
    cpp.append("          }\n")
    return cpp

