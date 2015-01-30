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
import subprocess
import hashlib
import shutil
import sys

class FunctionBucket:
  """A class that stores all the information necessary to write the ufl for a function (field or coefficient).
     Note that the class has limited ufl production because much of this is system dependent."""

  def __init__(self):
    """Define the expected members of the function class."""
    self.degree = None
    self.family = None
    self.name = None
    self.rank = None
    self.shape = None
    self.size = None
    self.enrichment_degree = None
    self.enrichment_family = None
    self.symbol = None
    self.symmetry = None
    self.type = None
    self.system = None
    self.functional = None
    self.cpp = None
    self.index = None
  
  def constant_ufl(self, suffix=""):
    """Returns a ufl string declaring a constant coefficient on a cell."""
    if self.rank == "Scalar":
      return self.symbol+suffix+" = Constant("+self.system.cell+")\n"
    elif self.rank == "Vector":
      return self.symbol+suffix+" = VectorConstant("+self.system.cell+")\n"
    elif self.rank == "Tensor":
      return self.symbol+suffix+" = TensorConstant("+self.system.cell+")\n"
    print self.rank
    print "Unknown rank."
    sys.exit(1)

  def element_ufl(self):
    ufl = []
    if self.enrichment_degree is not None and self.enrichment_family is not None:
      ufl = self.sub_element_ufl(index=0)
      ufl += self.sub_element_ufl(family=self.enrichment_family, degree=self.enrichment_degree, index=1)
      ufl.append(self.symbol+"_e = "+self.symbol+"_e0 + "+self.symbol+"_e1\n")
      ufl.append("\n")
    else:
      ufl = self.sub_element_ufl()
    return ufl

  def sub_element_ufl(self, family=None, degree=None, index=None):
    """Write an array of ufl strings describing the function (field or coefficient) element."""
    if family is None: family = self.family
    if degree is None: degree = self.degree
    if index is not None:
      index = `index`
    else:
      index = ""
    ufl = []
    ufl.append(declaration_comment("Element", self.type, self.name))
    ufl_line = self.symbol+"_e"+index+"= "
    if self.rank == "Scalar":
      ufl_line += "FiniteElement("
    elif self.rank == "Vector":
      # special case for vector elements - will need to be expanded
      if family in ["RT", "DRT"]:
        assert(self.size is None)
        ufl_line += "FiniteElement("
      else:
        ufl_line += "VectorElement("
    elif self.rank == "Tensor":
      ufl_line += "TensorElement("
    else:
      print self.rank
      print "Unknown rank."
      sys.exit(1)
    ufl_line += "\""+family +"\", " \
               +self.system.cell+", " \
               +`degree`
    if self.rank == "Vector":
      if self.size: ufl_line += ", dim="+`self.size`
    elif self.rank == "Tensor":
      if self.shape: ufl_line += ", shape=("+`self.shape[0]`+","+`self.shape[1]`+")"
      if self.symmetry: ufl_line += ", symmetry=True"
    ufl_line +=")\n"
    ufl.append(ufl_line)
    ufl.append("\n")
    return ufl

  def solvercoefficientspace_cpp(self, solvername, index=0, suffix=""):
    cpp = [] 
    if index == 0:
      cpp.append("        if (uflsymbol ==  \""+self.symbol+"\")\n")
    else:
      cpp.append("        else if (uflsymbol ==  \""+self.symbol+"\")\n")
    cpp.append("        {\n")
    cpp.append("          coefficientspace.reset(new "+self.system.name+solvername+"::CoefficientSpace_"+self.symbol+suffix+"(mesh));\n")
    cpp.append("        }\n")
    return cpp

  def namespace(self):
    return self.system.name+self.name

  def cppexpression_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the cpp expressions under this function."""
    cpp = []
    if index == 0:
      cpp.append("      if (functionname == \""+self.name+"\")\n")
    else:
      cpp.append("      else if (functionname ==  \""+self.name+"\")\n")
    cpp.append("      {\n")
    cpp.append("        if (expressiontype == \"initial_condition\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_cpp("initial_condition")
    cpp.append("        }\n")
    cpp.append("        else if (expressiontype == \"boundary_condition\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_cpp("boundary_condition")
    cpp.append("        }\n")
    cpp.append("        else if (expressiontype == \"value\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_cpp("value")
    cpp.append("        }\n")
    cpp.append("        else\n")
    cpp.append("        {\n")
    cpp.append("          tf_err(\"Unknown expressiontype in cpp_fetch_expression.\", \"Expression type: %s\", expressiontype.c_str());\n")
    cpp.append("        }\n")
    cpp.append("      }\n")
    return cpp

  def cppexpressiontype_cpp(self, basetype):
    """Write an array of cpp strings describing the namespace of the cpp expressions of a particular basetype under this function."""
    cpp = []

    expressions_found = 0
    for e in range(len(self.cpp)):
      if self.cpp[e].basetype==basetype:
        cpp += self.cpp[e].cppexpression_cpp(index=expressions_found)
        expressions_found += 1
    if expressions_found==0:
      cpp.append("          tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());\n")
    else:
      cpp.append("          else\n")
      cpp.append("          {\n")
      cpp.append("            tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());\n")
      cpp.append("          }\n")

    return cpp

  def cppexpression_init(self, index=0):
    """Write an array of cpp strings recasting an expression into the namespace of a cpp expression under this function."""
    cpp = []
    if index == 0:
      cpp.append("      if (functionname == \""+self.name+"\")\n")
    else:
      cpp.append("      else if (functionname ==  \""+self.name+"\")\n")
    cpp.append("      {\n")
    cpp.append("        if (expressiontype == \"initial_condition\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_init("initial_condition")
    cpp.append("        }\n")
    cpp.append("        else if (expressiontype == \"boundary_condition\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_init("boundary_condition")
    cpp.append("        }\n")
    cpp.append("        else if (expressiontype == \"value\")\n")
    cpp.append("        {\n")
    cpp += self.cppexpressiontype_init("value")
    cpp.append("        }\n")
    cpp.append("        else\n")
    cpp.append("        {\n")
    cpp.append("          tf_err(\"Unknown expressiontype in cpp_init_expression.\", \"Expression type: %s\", expressiontype.c_str());\n")
    cpp.append("        }\n")
    cpp.append("      }\n")
    return cpp

  def cppexpressiontype_init(self, basetype):
    """Write an array of cpp strings recasting an expression into the namespace of a cpp expression of a particular basetype under this function."""
    cpp = []

    expressions_found = 0
    for e in range(len(self.cpp)):
      if self.cpp[e].basetype==basetype:
        cpp += self.cpp[e].cppexpression_init(index=expressions_found)
        expressions_found += 1
    if expressions_found==0:
      cpp.append("          tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());\n")
    else:
      cpp.append("          else\n")
      cpp.append("          {\n")
      cpp.append("            tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());\n")
      cpp.append("          }\n")

    return cpp
