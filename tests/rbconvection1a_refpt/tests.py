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

stat = parser("rbconvection.stat")
det = parser("rbconvection.det")
steady = parser("rbconvection.steady")

def test_timestepcount():
  val = stat["timestep"]["value"][-1]
  assert abs(val - 35) < 1

def test_elapsedtime():
  val = stat["ElapsedTime"]["value"][-1]
  assert abs(val - 1900.0) < 5.e1

def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])
  assert abs(val - 42.865e-4) < 0.01

def test_nu():
  val = -1.0*(stat["Stokes"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

def test_extremumloc():
  val = (det["Stokes"]["Temperature"]["Array"][0:129/2,-1]).argmin()*1./128.
  assert abs(val - 0.2265625) < 0.01

def test_extremum():
  val = (det["Stokes"]["Temperature"]["Array"][0:129/2,-1]).min()
  assert abs(val - 0.4222) < 0.01

def test_steady():
  val = max(steady["Stokes"]["Velocity"]["L2NormSquared_change"][-1], \
            steady["Stokes"]["Pressure"]["Integral_change"][-1], \
            steady["Stokes"]["Dummy"]["Integral_change"][-1], \
            steady["Stokes"]["Temperature"]["BottomSurfaceIntegral_change"][-1], \
            steady["Stokes"]["Temperature"]["TopSurfaceIntegral_change"][-1])
  assert val < 1.e-5

def test_div():
  val = max(numpy.abs(stat["Divergence"]["Divergence"]["max"]).max(), \
            numpy.abs(stat["Divergence"]["Divergence"]["min"]).min())
  assert val < 1.e-6

def test_refpressure():
  val = det["Stokes"]["Pressure"]["Point"]
  assert numpy.all(val < 1.e-10)

