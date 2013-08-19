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


from buckettools.statfile import parser
import numpy
from math import sqrt

det = parser("madds1.det")

def test_P_0_1():
  val = det["Stokes"]["Pressure"]["Point"][0,-1]
  test = -1.2274872965e+00
  assert abs(val - test) < 0.001
  print "\tvalue=", val, "test=", test, " ",

def test_P_1_1():
  val = det["Stokes"]["Pressure"]["corner"][0,-1]
  test = 2.0861195420e-02
  assert abs(val - test) < 0.0001
  print "\tvalue=", val, "test=", test, " ",
