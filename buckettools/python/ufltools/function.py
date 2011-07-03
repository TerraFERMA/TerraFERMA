from ufltools.base import *

class Function:

  def __init__(self):
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
  
  def element_ufl(self):
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

