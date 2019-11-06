# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# (originally part of fluidity, modified for buckettools by Cian Wilson)

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


"""
Some mathematical routines
"""

import copy
import math
import unittest

import buckettools.diagnostics.debug as debug

try:
  import numpy
except ImportError:
  debug.deprint("Warning: Failed to import numpy module")
try:
  import scipy.linalg
except ImportError:
  debug.deprint("Warning: Failed to import scipy.linalg module")
try:
  import scipy.stats
except ImportError:
  debug.deprint("Warning: Failed to import scipy.stats module")

def Epsilon():
  """
  Return the smallest difference epsilon such that 1.0 + epsilon is
  distinguishable from 1.0
  """
  
  if not hasattr(Epsilon, "_epsilon"):
    # If we haven't yet calculated Epsilon, calculate it
    Epsilon._epsilon = 1.0
    for divide in [2.0, 1.75, 1.5, 1.25]:
      while True:
        newEpsilon = Epsilon._epsilon / divide
        if newEpsilon == Epsilon._epsilon:
          break
        elif 1.0 + newEpsilon > 1.0:
          Epsilon._epsilon = newEpsilon
        else:
          break
  
  return Epsilon._epsilon
  
def Nan():
  """
  Return floating point nan
  """
  
  return float("nan")
  
  
def AlmostEquals(x, y, tolerance = Epsilon()):
  """
  Return whether the supplied values are almost equal, to within the supplied
  tolerance
  """
  
  if abs(x) < tolerance:
    return abs(x - y) < tolerance
  else:
    return abs(x - y) / abs(x) < tolerance
    
