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
    self.symbol = None
    self.symmetry = None
    self.type = None
    self.system = None
    self.functional = None
    self.cpp = None
    self.functionals = None
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
    """Write an array of ufl strings describing the function (field or coefficient) element."""
    ufl = []
    ufl.append(declaration_comment("Element", self.type, self.name))
    ufl_line = self.symbol+"_e = "
    if self.rank == "Scalar":
      ufl_line += "FiniteElement("
    elif self.rank == "Vector":
      ufl_line += "VectorElement("
    elif self.rank == "Tensor":
      ufl_line += "TensorElement("
    else:
      print self.rank
      print "Unknown rank."
      sys.exit(1)
    ufl_line += "\""+self.family +"\", " \
               +self.system.cell+", " \
               +`self.degree`
    if self.rank == "Vector":
      if self.size: ufl_line += ", size="+`self.size`
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
    cpp.append("          dolfin::error(\"Unknown expressiontype in cpp_fetch_expression.\");\n")
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
      cpp.append("          dolfin::error(\"Unknown expressionname in cpp_fetch_expression.\");\n")
    else:
      cpp.append("          else\n")
      cpp.append("          {\n")
      cpp.append("            dolfin::error(\"Unknown expressionname in cpp_fetch_expression.\");\n")
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
    cpp.append("          dolfin::error(\"Unknown expressiontype in cpp_init_expression.\");\n")
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
      cpp.append("          dolfin::error(\"Unknown expressionname in cpp_fetch_expression.\");\n")
    else:
      cpp.append("          else\n")
      cpp.append("          {\n")
      cpp.append("            dolfin::error(\"Unknown expressionname in cpp_fetch_expression.\");\n")
      cpp.append("          }\n")

    return cpp
