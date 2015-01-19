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

class FunctionalBucket:
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl)."""

  def __init__(self):
    """Define the expected members of the functional class."""
    self.form = None
    self.system = None
    self.name = None
    self.symbol = None
    self.quadrature_degree = None
    self.quadrature_rule = None
  
  def ufl(self):
    """Write the functional to an array of ufl strings."""
    ufl = []
    for field in self.system.fields:
      ufl += field.element_ufl()
    ufl += self.system.element_ufl()
    ufl.append("\n")
    ufl += self.system.function_ufl()
    ufl.append("\n")
    ufl += self.system.iterate_ufl()
    ufl.append("\n")
    ufl += self.system.old_ufl()
    ufl.append("\n")
    for coeff in self.system.coeffs:
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
    for coeff in self.system.special_coeffs:
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
    for system in self.system.bucket.systems:
      if system.name == self.system.name: continue
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
    if self.system.bucket.parameters:
      ufl.append(comment("Global preamble"))
      ufl.append(self.system.bucket.parameters+"\n")
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
    ffc(self.namespace(), self.quadrature_rule, self.quadrature_degree)

  def namespace(self):
    name = self.system.name+self.name
    if self.function is not None: name += "_coefficient"
    return name

  def cpp(self, index=0):
    cpp = []
    if index == 0:
      cpp.append("        if (functionalname ==  \""+self.name+"\")\n")
    else:
      cpp.append("        else if (functionalname ==  \""+self.name+"\")\n")
    cpp.append("        {\n")
    cpp.append("          functional.reset(new "+self.namespace()+"::Form_"+self.symbol+"(mesh));\n")
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

  def repeated_uflsymbol_check(self):
    """Check for repeated ufl symbols."""
    stat = 0
    uflsymbols = self.system.bucket.list_globaluflsymbols()
    if self.symbol in uflsymbols:
      stat = 1
      print "ERROR functional ufl_symbol %s from functional %s::%s repeated in global ufl_symbols! Change one of its instances."%(self.symbol, self.system.name, self.name)
    repeated_auto_uflsymbols = set([s for s in uflsymbols for a in uflsymbol_suffixes() if s+a == self.symbol and a != ''])
    if len(repeated_auto_uflsymbols) > 0: stat = 1
    for s in repeated_auto_uflsymbols: print "ERROR functional ufl_symbol %s from functional %s::%s matches ufl_symbol generated from global ufl_symbol %s! Change functional ufl_symbol %s to avoid reserved endings."%(self.symbol, self.system.name, self.name, s, self.symbol)
    return stat

