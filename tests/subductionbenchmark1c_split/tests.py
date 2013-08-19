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

det = parser("subduction.det")

def test_T_11_11():
  val = det["Temperature"]["Temperature"]["SlabPoint"][0,-1]-273.
  test=388.5
  assert abs(val - 388.5) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_slab():
  val = sqrt(sum((det["Temperature"]["Temperature"]["Slab"][:,-1]-273.)**2)/36.)
  test= 504.0
  assert abs(val - test) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_wedge():
  val = sqrt(sum((det["Temperature"]["Temperature"]["Wedge"][:,-1]-273.)**2)/78.)
  test=853.5
  assert abs(val - 853.5) < 1.0
  print "\tvalue=", val, "test=", test, " ",

