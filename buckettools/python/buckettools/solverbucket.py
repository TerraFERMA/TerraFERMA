# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

from buckettools.base import *
import sys
import subprocess
import hashlib
import shutil

class SolverBucket:
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
    self.quadrature_degree = None
    self.quadrature_rule = None

  def ufl(self):
    """Write the system of forms to an array."""
    ufl = []
    ufl += self.system.functions_ufl()
    
    if self.system.bucket.parameters:
      ufl.append(comment("Global preamble"))
      ufl.append(self.system.bucket.parameters+"\n")
    
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
    if suffix: filename += suffix
    filehandle = file(filename, 'w')
    filehandle.writelines(ufl)
    filehandle.close()
    
  def write_ufc(self):
    """Write the system of forms to a ufl file."""
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

      qri  = headertext.find("quadrature_rule")
      qrin = headertext.find("\n", qri, -1)
      if (qri != -1) and (qrin != -1):
        qr_old = headertext[qri:qrin].split(" ")[-1]
        qr_new = "'"+self.quadrature_rule+"'" 
        rebuild = rebuild or qr_new != qr_old
      else: # if we don't know what the old quadrature rule was we always want to (re)build
        rebuild = True
      
    else:   # if we don't have a header file we always want to (re)build
      rebuild = True

    if rebuild:
      # files and/or quadrature degree have changed
      shutil.copy(filename+".temp", filename)
      command = ["ffc", "-l", "dolfin", "-O", "-r", "quadrature"]
      if self.quadrature_degree:
        command += ["-fquadrature_degree="+`self.quadrature_degree`]
      command += ["-q", self.quadrature_rule]
      command += [filename]
      try:
        subprocess.check_call(command)
      except:
        print "ERROR while calling ffc on file ", filename
        sys.exit(1)
    
  def functionspace_cpp_no_if(self):
    cpp = []
    cpp.append("      if (periodicmap && facetdomains)\n")
    cpp.append("      {\n")
    cpp.append("        functionspace.reset( new "+self.namespace()+"::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );\n")
    cpp.append("      }\n")
    cpp.append("      else\n")
    cpp.append("      {\n")
    cpp.append("        functionspace.reset( new "+self.namespace()+"::FunctionSpace(mesh) );\n")
    cpp.append("      }\n")
    return cpp

  def functionspace_cpp(self, index=0):
    cpp = [] 
    if index == 0:
      cpp.append("      if (solvername ==  \""+self.name+"\")\n")
    else:
      cpp.append("      else if (solvername ==  \""+self.name+"\")\n")
    cpp.append("      {\n")
    cpp.append("        if (periodicmap && facetdomains)\n")
    cpp.append("        {\n")
    cpp.append("          functionspace.reset(new "+self.namespace()+"::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );\n") 
    cpp.append("        }\n")
    cpp.append("        else\n")
    cpp.append("        {\n")
    cpp.append("          functionspace.reset(new "+self.namespace()+"::FunctionSpace(mesh));\n")
    cpp.append("        }\n")
    cpp.append("      }\n")
    return cpp

  def form_cpp(self):
    cpp = []
    for f in range(len(self.forms)):
      if f == 0:
        cpp.append("          if (formname == \""+self.form_names[f]+"\")\n")
      else:
        cpp.append("          else if (formname == \""+self.form_names[f]+"\")\n")
      cpp.append("          {\n")
      if self.form_ranks[f]==0:
        cpp.append("            form.reset(new "+self.system.name+self.name+"::Form_"+self.form_symbols[f]+"(functionspace));\n")
      elif self.form_ranks[f]==1:
        cpp.append("            form.reset(new "+self.system.name+self.name+"::Form_"+self.form_symbols[f]+"(functionspace, functionspace));\n")
      else:
        print "Unknwon form rank."
        sys.exit(1)
      cpp.append("          }\n")
    return cpp
      

