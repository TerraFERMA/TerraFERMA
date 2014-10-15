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

class CppExpressionBucket:
  """A class that stores all the information necessary to write the cpp for a user defined expression."""

  def __init__(self):
    self.name     = None
    self.members  = None
    self.initfunc = None
    self.evalfunc = None
    self.function = None
    self.rank     = None
    self.basetype = None
    self.nametype = None
    self.include  = None
  
  def namespace(self):
    return self.function.system.name+self.function.name+self.name+self.nametype+"Expression"

  def cpp(self):
    """Write the cpp expression to an array of cpp header strings."""
    cpp = []

    cpp.append("#ifndef __"+self.namespace().upper()+"_H\n")
    cpp.append("#define __"+self.namespace().upper()+"_H\n")
    cpp.append("\n")
    cpp.append("#include \"BoostTypes.h\"\n")
    cpp.append("#include \"Bucket.h\"\n")
    cpp.append("#include \"SystemBucket.h\"\n")
    cpp.append("#include <dolfin.h>\n")
    if self.include:
      for line in self.include.split("\n"):
        cpp.append(line+"\n")
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
      cpp.append("    "+self.namespace()+"(const Bucket *bucket, const SystemBucket *system, const double_ptr time) : dolfin::Expression(), bucket_(bucket), system_(system), time_(time), initialized_(false)\n")
    elif self.rank == "Vector":
      cpp.append("    "+self.namespace()+"(const std::size_t &dim, const Bucket *bucket, const SystemBucket *system, const double_ptr time) : dolfin::Expression(dim), bucket_(bucket), system_(system), time_(time), initialized_(false)\n")
    elif self.rank == "Tensor":
      cpp.append("    "+self.namespace()+"(const std::vector<std::size_t> &value_shape, const Bucket *bucket, const SystemBucket *system, const double_ptr time) : dolfin::Expression(value_shape), bucket_(bucket), system_(system), time_(time), initialized_(false)\n")
    else:
      print self.rank
      print "Unknown rank."
      sys.exit(1)
    cpp.append("    {\n")
    cpp.append("                                                                     // do nothing\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const\n")
    cpp.append("    {\n")
    cpp.append("      \n")
    cpp.append("      tf_err(\"Buckettools C++ expressions must be called using the eval(values, x, cell) interface.\", \"Cannot use eval(values, x) interface.\");\n")
    cpp.append("      \n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const\n")
    cpp.append("    {\n")
    cpp.append("      \n")
    for line in self.evalfunc.split("\n"):
      cpp.append("      "+line+"\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    void init()\n")
    cpp.append("    {\n")
    cpp.append("      if (!initialized_)\n")
    cpp.append("      {\n")
    cpp.append("        initialized_ = true;\n")
    cpp.append("        \n")
    for line in self.initfunc.split("\n"):
      cpp.append("        "+line+"\n")
    cpp.append("      }\n")
    cpp.append("    }\n")
    cpp.append("  \n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  // Private functions\n")
    cpp.append("  //*****************************************************************|***********************************************************//\n")
    cpp.append("  \n")
    cpp.append("  private:                                                           // only available to this class\n")
    cpp.append("    \n")
    cpp.append("    const Bucket* bucket() const\n")
    cpp.append("    {\n")
    cpp.append("      return bucket_;\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    const SystemBucket* system() const\n")
    cpp.append("    {\n")
    cpp.append("      return system_;\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    const double_ptr time() const\n")
    cpp.append("    {\n")
    cpp.append("      return time_;\n")
    cpp.append("    }\n")
    cpp.append("    \n")
    cpp.append("    const Bucket *bucket_;\n")
    cpp.append("    \n")
    cpp.append("    const SystemBucket *system_;\n")
    cpp.append("    \n")
    cpp.append("    const double_ptr time_;\n")
    cpp.append("    \n")
    cpp.append("    bool initialized_;\n")
    cpp.append("    \n")
    for line in self.members.split("\n"):
      cpp.append("    "+line+"\n")
    cpp.append("  \n")
    cpp.append("  };\n")
    cpp.append("  \n")
    cpp.append("}\n")

    cpp.append("\n")
    cpp.append("#endif\n")
    cpp.append("\n")

    return cpp

  def write_cppheader(self, suffix=None):
    """Write the cpp expression to a cpp header file."""
    cpp = self.cpp()

    filename   = self.namespace()+".h"
    if suffix: filename += suffix 
    filehandle = file(filename, 'w')
    filehandle.writelines(cpp)
    filehandle.close()

  def write_cppheader_md5(self):
    """Write the cpp expression to a cpp header file (with md5 checksum)."""
    self.write_cppheader(suffix=".temp")

    filename = self.namespace()+".h"

    try:
      checksum = hashlib.md5(open(filename).read()).hexdigest()
    except:
      checksum = None

    if checksum != hashlib.md5(open(filename+".temp").read()).hexdigest():
      # files have changed
      shutil.copy(filename+".temp", filename)

  def cppexpression_cpp(self, index=0):
    """Write an array of cpp strings describing the namespace of the cpp expression."""
    cpp = []
    
    if index == 0:
      cpp.append("          if (expressionname == \""+self.name+"\")\n")
    else:
      cpp.append("          else if (expressionname ==  \""+self.name+"\")\n")
    cpp.append("          {\n")
    if self.rank == "Scalar":
      cpp.append("            expression.reset(new "+self.namespace()+"(bucket, system, time));\n")
    elif self.rank == "Vector":
      cpp.append("            expression.reset(new "+self.namespace()+"(size, bucket, system, time));\n")
    elif self.rank == "Tensor":
      cpp.append("            expression.reset(new "+self.namespace()+"(shape, bucket, system, time));\n")
    else:
      print self.rank
      print "Unknown rank."
      sys.exit(1)
    cpp.append("          }\n")
    
    return cpp

  def cppexpression_init(self, index=0):
    """Write an array of cpp strings recasting an expression into the namespace of the cpp expression."""
    cpp = []
    
    if index == 0:
      cpp.append("          if (expressionname == \""+self.name+"\")\n")
    else:
      cpp.append("          else if (expressionname ==  \""+self.name+"\")\n")
    cpp.append("          {\n")
    cpp.append("            (*std::dynamic_pointer_cast< "+self.namespace()+" >(expression)).init();\n")
    cpp.append("          }\n")
    
    return cpp

