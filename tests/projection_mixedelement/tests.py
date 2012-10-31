
from buckettools.statfile import parser
import numpy

stat = parser("projection.stat")

def test_picard_field1_min():
  val = stat["PicardProjection"]["Field1"]["min"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_field1_max():
  val = stat["PicardProjection"]["Field1"]["max"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_field1_int():
  val = stat["PicardProjection"]["Field1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_source1_int():
  val = stat["PicardProjection"]["Source1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_picard_field2_min_0():
  val = stat["PicardProjection"]["Field2"]["max"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_field2_max_0():
  val = stat["PicardProjection"]["Field2"]["min"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_field2_int_0():
  val = stat["PicardProjection"]["Field2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_source2_int_0():
  val = stat["PicardProjection"]["Source2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_picard_field2_min_1():
  val = stat["PicardProjection"]["Field2"]["max"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_field2_max_1():
  val = stat["PicardProjection"]["Field2"]["min"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_field2_int_1():
  val = stat["PicardProjection"]["Field2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_source2_int_1():
  val = stat["PicardProjection"]["Source2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_picard_dummy_max():
  val = stat["PicardProjection"]["Dummy"]["max"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_picard_dummy_min():
  val = stat["PicardProjection"]["Dummy"]["min"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_snes_field1_min():
  val = stat["SNESProjection"]["Field1"]["min"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_field1_max():
  val = stat["SNESProjection"]["Field1"]["max"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_field1_int():
  val = stat["SNESProjection"]["Field1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_source1_int():
  val = stat["SNESProjection"]["Source1"]["Integral"][-1]
  assert abs(val - 1000.0) < 1.e-7

def test_snes_field2_min_0():
  val = stat["SNESProjection"]["Field2"]["max"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_field2_max_0():
  val = stat["SNESProjection"]["Field2"]["min"][0][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_field2_int_0():
  val = stat["SNESProjection"]["Field2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_source2_int_0():
  val = stat["SNESProjection"]["Source2"]["Integral0"][-1]
  assert abs(val - 2000.0) < 1.e-7

def test_snes_field2_min_1():
  val = stat["SNESProjection"]["Field2"]["max"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_field2_max_1():
  val = stat["SNESProjection"]["Field2"]["min"][1][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_field2_int_1():
  val = stat["SNESProjection"]["Field2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_source2_int_1():
  val = stat["SNESProjection"]["Source2"]["Integral1"][-1]
  assert abs(val - 3000.0) < 1.e-7

def test_snes_dummy_max():
  val = stat["SNESProjection"]["Dummy"]["max"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

def test_snes_dummy_min():
  val = stat["SNESProjection"]["Dummy"]["min"]
  assert numpy.all(abs(val - 6.0) < 1.e-7)

