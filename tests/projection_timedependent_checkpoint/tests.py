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
import libspud

stat = parser("projection_checkpoint.stat")

def test_picard_field1_min():
  val = stat["PicardProjection"]["Field1"]["min"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_field1_max():
  val = stat["PicardProjection"]["Field1"]["max"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_field1_int():
  val = stat["PicardProjection"]["Field1"]["Integral"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_source1_int():
  val = stat["PicardProjection"]["Source1"]["Integral"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_field1_oldint():
  val = stat["PicardProjection"]["Field1"]["OldIntegral"]
  assert numpy.all(val == [500.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_picard_source1_oldint():
  val = stat["PicardProjection"]["Source1"]["OldIntegral"]
  assert numpy.all(val == [500.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_picard_field2_min_0():
  val = stat["PicardProjection"]["Field2"]["max"][0]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_field2_max_0():
  val = stat["PicardProjection"]["Field2"]["min"][0]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_field2_int_0():
  val = stat["PicardProjection"]["Field2"]["Integral0"]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_source2_int_0():
  val = stat["PicardProjection"]["Source2"]["Integral0"]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_field2_oldint_0():
  val = stat["PicardProjection"]["Field2"]["OldIntegral0"]
  assert numpy.all(val == [1000.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_picard_source2_oldint_0():
  val = stat["PicardProjection"]["Source2"]["OldIntegral0"]
  assert numpy.all(val == [1000.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_picard_field2_min_1():
  val = stat["PicardProjection"]["Field2"]["max"][1]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_field2_max_1():
  val = stat["PicardProjection"]["Field2"]["min"][1]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_field2_int_1():
  val = stat["PicardProjection"]["Field2"]["Integral1"]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_source2_int_1():
  val = stat["PicardProjection"]["Source2"]["Integral1"]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_field2_oldint_1():
  val = stat["PicardProjection"]["Field2"]["OldIntegral1"]
  assert numpy.all(val == [1500.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_picard_source2_oldint_1():
  val = stat["PicardProjection"]["Source2"]["OldIntegral1"]
  assert numpy.all(val == [1500.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_snes_field1_min():
  val = stat["SNESProjection"]["Field1"]["min"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_max():
  val = stat["SNESProjection"]["Field1"]["max"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_int():
  val = stat["SNESProjection"]["Field1"]["Integral"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_source1_int():
  val = stat["SNESProjection"]["Source1"]["Integral"]
  assert numpy.all(val == [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_oldint():
  val = stat["SNESProjection"]["Field1"]["OldIntegral"]
  assert numpy.all(val == [500.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_snes_source1_oldint():
  val = stat["SNESProjection"]["Source1"]["OldIntegral"]
  assert numpy.all(val == [500.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_snes_field2_min_0():
  val = stat["SNESProjection"]["Field2"]["max"][0]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_max_0():
  val = stat["SNESProjection"]["Field2"]["min"][0]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_int_0():
  val = stat["SNESProjection"]["Field2"]["Integral0"]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_source2_int_0():
  val = stat["SNESProjection"]["Source2"]["Integral0"]
  assert numpy.all(val == [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_oldint_0():
  val = stat["SNESProjection"]["Field2"]["OldIntegral0"]
  assert numpy.all(val == [1000.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_snes_source2_oldint_0():
  val = stat["SNESProjection"]["Source2"]["OldIntegral0"]
  assert numpy.all(val == [1000.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_snes_field2_min_1():
  val = stat["SNESProjection"]["Field2"]["max"][1]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_max_1():
  val = stat["SNESProjection"]["Field2"]["min"][1]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_int_1():
  val = stat["SNESProjection"]["Field2"]["Integral1"]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_source2_int_1():
  val = stat["SNESProjection"]["Source2"]["Integral1"]
  assert numpy.all(val == [1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_oldint_1():
  val = stat["SNESProjection"]["Field2"]["OldIntegral1"]
  assert numpy.all(val == [1500.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_snes_source2_oldint_1():
  val = stat["SNESProjection"]["Source2"]["OldIntegral1"]
  assert numpy.all(val == [1500.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_checkpoint_0():
  libspud.load_options("projection_checkpoint_checkpoint_0.tfml")
  val = libspud.get_option("/timestepping/current_time")
  libspud.clear_options()
  assert val==10.0

