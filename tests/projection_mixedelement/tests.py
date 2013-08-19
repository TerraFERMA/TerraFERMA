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
import buckettools.vtktools as vtktools
import numpy

vtu = vtktools.vtu("projection000010.vtu")
stat = parser("projection.stat")
locations = numpy.array([[ 0.5 ,  0.  ,  0.  ],
 [ 0.25,  0.25,  0.  ],
 [ 0.  ,  0.  ,  0.  ],
 [ 0.  ,  0.5 ,  0.  ],
 [ 0.75,  0.25,  0.  ],
 [ 1.  ,  0.  ,  0.  ],
 [ 1.  ,  0.5 ,  0.  ],
 [ 0.5 ,  0.5 ,  0.  ],
 [ 0.75,  0.75,  0.  ],
 [ 0.25,  0.75,  0.  ],
 [ 0.  ,  1.  ,  0.  ],
 [ 1.  ,  1.  ,  0.  ],
 [ 0.5 ,  1.  ,  0.  ]]) # the p2 locations


def test_timestep_count():
  val = stat["timestep"]["value"][-1]
  assert val==10

def test_picard_field1_min():
  val = stat["PicardProjection"]["Field1"]["min"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_field1_max():
  val = stat["PicardProjection"]["Field1"]["max"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_field1_vtu():
  val = vtu.GetField("PicardProjection::Field1")
  assert numpy.all(abs(val-1000.0) < 1.e-9)

def test_picard_resfield1_vtu():
  val = vtu.GetField("PicardProjection::ResidualField1")
  assert numpy.all(abs(val) < 1.e-9)

def test_picard_field1_int():
  val = stat["PicardProjection"]["Field1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_source1_int():
  val = stat["PicardProjection"]["Source1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_source1_vtu():
  val = vtu.GetField("PicardProjection::Source1")
  assert numpy.all(abs(val-1000.0) < 1.e-9)

def test_picard_field2_min_0():
  val = stat["PicardProjection"]["Field2"]["max"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_field2_max_0():
  val = stat["PicardProjection"]["Field2"]["min"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_field2_vtu_0():
  val = vtu.GetField("PicardProjection::Field2")[:,0]
  assert numpy.all(abs(val-2000.0) < 1.e-9)

def test_picard_resfield2_vtu_0():
  val = vtu.GetField("PicardProjection::ResidualField2")[:,0]
  assert numpy.all(abs(val) < 1.e-9)

def test_picard_field2_int_0():
  val = stat["PicardProjection"]["Field2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_source2_int_0():
  val = stat["PicardProjection"]["Source2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_source2_vtu_0():
  val = vtu.GetField("PicardProjection::Source2")[:,0]
  assert numpy.all(abs(val-2000.0) < 1.e-9)

def test_picard_field2_min_1():
  val = stat["PicardProjection"]["Field2"]["max"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_field2_max_1():
  val = stat["PicardProjection"]["Field2"]["min"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_field2_vtu_1():
  val = vtu.GetField("PicardProjection::Field2")[:,1]
  assert numpy.all(abs(val-3000.0) < 1.e-9)

def test_picard_resfield2_vtu_1():
  val = vtu.GetField("PicardProjection::ResidualField2")[:,1]
  assert numpy.all(abs(val) < 1.e-9)

def test_picard_field2_int_1():
  val = stat["PicardProjection"]["Field2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_source2_int_1():
  val = stat["PicardProjection"]["Source2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_source2_vtu_1():
  val = vtu.GetField("PicardProjection::Source2")[:,1]
  assert numpy.all(abs(val-3000.0) < 1.e-9)

def test_picard_dummy_max():
  val = stat["PicardProjection"]["Dummy"]["max"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_picard_dummy_min():
  val = stat["PicardProjection"]["Dummy"]["min"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_picard_dummy_vtu():
  val = vtu.GetField("PicardProjection::Dummy")
  assert numpy.all(abs(val-6.0) < 1.e-9)

def test_snes_field1_min():
  val = stat["SNESProjection"]["Field1"]["min"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_field1_max():
  val = stat["SNESProjection"]["Field1"]["max"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_field1_vtu_1():
  val = vtu.GetField("SNESProjection::Field1")
  assert numpy.all(abs(val-1000.0) < 1.e-9)

def test_snes_resfield1_vtu():
  val = vtu.GetField("SNESProjection::ResidualField1")
  assert numpy.all(abs(val) < 1.e-9)

def test_snes_field1_int():
  val = stat["SNESProjection"]["Field1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_source1_int():
  val = stat["SNESProjection"]["Source1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_source1_vtu():
  val = vtu.GetField("SNESProjection::Source1")[:,0]
  assert numpy.all(abs(val-1000.0) < 1.e-9)

def test_snes_field2_min_0():
  val = stat["SNESProjection"]["Field2"]["max"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_field2_max_0():
  val = stat["SNESProjection"]["Field2"]["min"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_field2_vtu_0():
  val = vtu.GetField("SNESProjection::Field2")[:,0]
  assert numpy.all(abs(val-2000.0) < 1.e-9)

def test_snes_resfield2_vtu_0():
  val = vtu.GetField("SNESProjection::ResidualField2")[:,0]
  assert numpy.all(abs(val) < 1.e-9)

def test_snes_field2_int_0():
  val = stat["SNESProjection"]["Field2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_source2_int_0():
  val = stat["SNESProjection"]["Source2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_source2_vtu_0():
  val = vtu.GetField("SNESProjection::Source2")[:,0]
  assert numpy.all(abs(val-2000.0) < 1.e-9)

def test_snes_field2_min_1():
  val = stat["SNESProjection"]["Field2"]["max"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_field2_max_1():
  val = stat["SNESProjection"]["Field2"]["min"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_field2_vtu_1():
  val = vtu.GetField("SNESProjection::Field2")[:,1]
  assert numpy.all(abs(val-3000.0) < 1.e-9)

def test_snes_resfield2_vtu_1():
  val = vtu.GetField("SNESProjection::ResidualField2")[:,1]
  assert numpy.all(abs(val) < 1.e-9)

def test_snes_field2_int_1():
  val = stat["SNESProjection"]["Field2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_source2_int_1():
  val = stat["SNESProjection"]["Source2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_source2_vtu_1():
  val = vtu.GetField("SNESProjection::Source2")[:,1]
  assert numpy.all(abs(val-3000.0) < 1.e-9)

def test_snes_dummy_max():
  val = stat["SNESProjection"]["Dummy"]["max"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_snes_dummy_min():
  val = stat["SNESProjection"]["Dummy"]["min"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_snes_dummy_vtu():
  val = vtu.GetField("SNESProjection::Dummy")
  assert numpy.all(abs(val-6.0) < 1.e-9)

def test_locations_vtu():
  val = vtu.GetLocations()
  assert numpy.all(val == locations)

