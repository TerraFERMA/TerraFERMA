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

import libspud
import buckettools.cppexpressionbucket
import sys

class SpudCppExpressionBucket(buckettools.cppexpressionbucket.CppExpressionBucket):
  """A class that stores all the information necessary to write the cpp for a user defined expression 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, name, function):
    """Fill a cpp expression class with data describing that expression using libspud, the given optionpath and the function its based on."""
    self.name     = name
    self.members  = libspud.get_option(optionpath+"/cpp/members")
    self.initfunc = libspud.get_option(optionpath+"/cpp/initialization")
    self.evalfunc = libspud.get_option(optionpath+"/cpp/eval")
    self.function = function
    
    if libspud.have_option(optionpath+"/cpp/include"):
      self.include = libspud.get_option(optionpath+"/cpp/include")

    self.basetype     = libspud.get_option(optionpath+"/type")
    if self.basetype=="initial_condition":
      self.nametype = "IC"
    elif self.basetype=="value":
      self.nametype = "Value"
    elif self.basetype=="boundary_condition":
      self.nametype = "BC"
    else:
      print self.basetype
      print "Unknown type."
      sys.exit(1)

    rank          = libspud.get_option(optionpath+"/cpp/rank")
    if rank=="0": # for consistency with other parts of the code convert this to a human parseable string
      self.rank = "Scalar"
    elif rank=="1":
      self.rank = "Vector"
    elif rank=="2":
      self.rank = "Tensor"
    else:
      print rank
      print "Unknown rank."
      sys.exit(1)
      

