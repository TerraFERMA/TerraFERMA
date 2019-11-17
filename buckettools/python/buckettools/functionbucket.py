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
import os

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
    self.quadrature_rule = None
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
      return self.symbol+suffix+" = Constant("+self.system.cell+")"+os.linesep
    elif self.rank == "Vector":
      return self.symbol+suffix+" = VectorConstant("+self.system.cell+")"+os.linesep
    elif self.rank == "Tensor":
      return self.symbol+suffix+" = TensorConstant("+self.system.cell+")"+os.linesep
    print(self.rank)
    print("Unknown rank.")
    sys.exit(1)

  def element_ufl(self):
    if self.enrichment_degree is not None and self.enrichment_family is not None:
      ufl = self.scalar_element_ufl(index=0)
      ufl += self.scalar_element_ufl(family=self.enrichment_family, degree=self.enrichment_degree, index=1)
      ufl.append(self.symbol+"_es = "+self.symbol+"_es0 + "+self.symbol+"_es1"+os.linesep)
      ufl.append(os.linesep)
    else:
      ufl = self.scalar_element_ufl()
    
    ufl_line = self.symbol+"_e = "
    if self.rank == "Scalar" or (self.rank == "Vector" and self.family in ["RT", "DRT", "BDM", "N1curl", "N2curl"]):
      assert(self.size is None)
      ufl_line += self.symbol+"_es"
    else:
      if self.rank == "Vector":
        ufl_line += "VectorElement("+self.symbol+"_es"
        if self.size: ufl_line += ", dim="+repr(self.size)
      elif self.rank == "Tensor":
        ufl_line += "TensorElement("+self.symbol+"_es"
        if self.shape: ufl_line += ", shape=("+repr(self.shape[0])+","+repr(self.shape[1])+")"
        if self.symmetry: ufl_line += ", symmetry=True"
      else:
        print(self.rank)
        print("Unknown rank.")
        sys.exit(1)
      ufl_line +=")"+os.linesep
    ufl.append(ufl_line)
    ufl.append(os.linesep)
    return ufl

  def scalar_element_ufl(self, family=None, degree=None, index=None):
    """Write an array of ufl strings describing the function (field or coefficient) element."""
    if family is None: family = self.family
    if degree is None: degree = self.degree
    if index is not None:
      index = repr(index)
    else:
      index = ""
    ufl = []
    ufl.append(declaration_comment("Element", self.type, self.name))
    ufl_line = self.symbol+"_es"+index+" = "
    ufl_line += "FiniteElement("
    ufl_line += "\""+family +"\", " \
               +self.system.cell+", " \
               +repr(degree)
    if self.quadrature_rule is not None: ufl_line += ", quad_scheme="+"\""+self.quadrature_rule+"\""
    ufl_line +=")"+os.linesep
    ufl.append(ufl_line)
    ufl.append(os.linesep)
    return ufl

  def solvercoefficientspace_cpp(self, solvername, index=0, suffix=""):
    cpp = [] 
    if index == 0:
      cpp.append("        if (uflsymbol ==  \""+self.symbol+"\")"+os.linesep)
    else:
      cpp.append("        else if (uflsymbol ==  \""+self.symbol+"\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp.append("          coefficientspace.reset(new "+self.system.name+solvername+"::CoefficientSpace_"+self.symbol+suffix+"(mesh));"+os.linesep)
    cpp.append("        }"+os.linesep)
    return cpp

  def cppexpression_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the cpp expressions under this function."""
    cpp = []
    if index == 0:
      cpp.append("      if (functionname == \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("      else if (functionname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("      {"+os.linesep)
    cpp.append("        if (expressiontype == \"initial_condition\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_cpp("initial_condition")
    cpp.append("        }"+os.linesep)
    cpp.append("        else if (expressiontype == \"boundary_condition\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_cpp("boundary_condition")
    cpp.append("        }"+os.linesep)
    cpp.append("        else if (expressiontype == \"value\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_cpp("value")
    cpp.append("        }"+os.linesep)
    cpp.append("        else"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp.append("          tf_err(\"Unknown expressiontype in cpp_fetch_expression.\", \"Expression type: %s\", expressiontype.c_str());"+os.linesep)
    cpp.append("        }"+os.linesep)
    cpp.append("      }"+os.linesep)
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
      cpp.append("          tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());"+os.linesep)
    else:
      cpp.append("          else"+os.linesep)
      cpp.append("          {"+os.linesep)
      cpp.append("            tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("          }"+os.linesep)

    return cpp

  def cppexpression_init(self, index=0):
    """Write an array of cpp strings recasting an expression into the namespace of a cpp expression under this function."""
    cpp = []
    if index == 0:
      cpp.append("      if (functionname == \""+self.name+"\")"+os.linesep)
    else:
      cpp.append("      else if (functionname ==  \""+self.name+"\")"+os.linesep)
    cpp.append("      {"+os.linesep)
    cpp.append("        if (expressiontype == \"initial_condition\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_init("initial_condition")
    cpp.append("        }"+os.linesep)
    cpp.append("        else if (expressiontype == \"boundary_condition\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_init("boundary_condition")
    cpp.append("        }"+os.linesep)
    cpp.append("        else if (expressiontype == \"value\")"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp += self.cppexpressiontype_init("value")
    cpp.append("        }"+os.linesep)
    cpp.append("        else"+os.linesep)
    cpp.append("        {"+os.linesep)
    cpp.append("          tf_err(\"Unknown expressiontype in cpp_init_expression.\", \"Expression type: %s\", expressiontype.c_str());"+os.linesep)
    cpp.append("        }"+os.linesep)
    cpp.append("      }"+os.linesep)
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
      cpp.append("          tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());"+os.linesep)
    else:
      cpp.append("          else"+os.linesep)
      cpp.append("          {"+os.linesep)
      cpp.append("            tf_err(\"Unknown expressionname in cpp_fetch_expression.\", \"Expression name: %s\", functionname.c_str());"+os.linesep)
      cpp.append("          }"+os.linesep)

    return cpp
