from ufltools.base import *
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

  def cppexpression(self):
    """Write the cpp expression to an array of cpp header strings."""
    assert(self.cpp)
    
    cpp = []
    cpp.append("#ifndef __"+self.namespace().upper()+"_EXPRESSION_H\n")
    cpp.append("#define __"+self.namespace().upper()+"_EXPRESSION_H\n")
    cpp.append("\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include \"Bucket.h\"\n")
    cpp.append("#include <dolfin.h>\n")
    cpp.append("\n")

    cpp.append("namespace buckettools\n")
    cpp.append("{\n")
    cpp.append("  //*****************************************************************|************************************************************//\n")
    cpp.append("  // "+self.namespace()+" class:\n")
    cpp.append("  //\n")
    cpp.append("  // The "+self.namespace()+" class describes a derived dolfin Expression class that overloads\n")
    cpp.append("  // the eval function using a user defined data.\n")
    cpp.append("  //*****************************************************************|************************************************************//\n")
    cpp.append("  class "+self.namespace()+" : public dolfin::Expression\n")
    cpp.append("  {\n")
    cpp.append("  \n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  // Publicly available functions\n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  \n")
    cpp.append("  public:                                                            // available to everyone\n")
    cpp.append("  \n")
    if self.rank == "Scalar":
      cpp.append("    "+self.namespace()+"(const Bucket *bucket) : dolfin::Expression(), bucket(bucket)\n")
    elif self.rank == "Vector":
      cpp.append("    "+self.namespace()+"(const uint &dim, const Bucket *bucket) : dolfin::Expression(dim), bucket(bucket)\n")
    elif self.rank == "Tensor":
      cpp.append("    "+self.namespace()+"(const std::vector<uint> &value_shape, const Bucket *bucket) : dolfin::Expression(value_shape), bucket(bucket)\n")
    else:
      print self.rank
      print "Unknown rank."
      sys.exit(1)
    cpp.append("    {\n")
    for line in self.cpp["constructor"].split("\n"):
      cpp.append("      "+line+"\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell)\n")
    cpp.append("    {\n")
    for line in self.cpp["eval"].split("\n"):
      cpp.append("      "+line+"\n")
    cpp.append("    }\n")
    cpp.append("  \n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  // Private functions\n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  \n")
    cpp.append("  private:                                                           // only available to this class\n")
    cpp.append("    \n")
    cpp.append("    const Bucket *bucket;\n")
    cpp.append("    \n")
    for line in self.cpp["members"].split("\n"):
      cpp.append("    "+line+"\n")
    cpp.append("  \n")
    cpp.append("  };\n")
    cpp.append("  \n")
    cpp.append("}\n")

    cpp.append("\n")
    cpp.append("#endif\n")
    cpp.append("\n")

    return cpp

  def write_cppexpression(self, suffix=None):
    """Write the cpp expression to a cpp header file."""
    cpp = self.cppexpression()

    filename   = self.namespace()+".h"
    if suffix: filename += suffix 
    filehandle = file(filename, 'w')
    filehandle.writelines(cpp)
    filehandle.close()

  def write_cppexpressionheader(self):
    """Write the cpp expression to a cpp header file (with md5 checksum)."""
    self.write_cppexpression(suffix=".temp")

    filename = self.namespace()+".h"

    try:
      checksum = hashlib.md5(open(filename).read()).hexdigest()
    except:
      checksum = None

    if checksum != hashlib.md5(open(filename+".temp").read()).hexdigest():
      # files have changed
      shutil.copy(filename+".temp", filename)

