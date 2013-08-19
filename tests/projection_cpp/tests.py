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

stat = parser("projection.stat")

def test_picard_field1_min():
  val = stat["PicardProjection"]["Field1"]["min"]
  assert numpy.all(val == [0.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_field1_max():
  val = stat["PicardProjection"]["Field1"]["max"]
  assert numpy.all(val == [2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_field1_int():
  val = stat["PicardProjection"]["Field1"]["Integral"]
  assert numpy.all(abs(val - [1.33333333, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0]) < 1.e-6)

def test_picard_source1_min():
  val = stat["PicardProjection"]["Source1"]["min"]
  assert numpy.all(val == [2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_source1_max():
  val = stat["PicardProjection"]["Source1"]["max"]
  assert numpy.all(val == [2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_source1_int():
  val = stat["PicardProjection"]["Source1"]["Integral"]
  assert numpy.all(val == [2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_pysource1_int():
  val = stat["PicardProjection"]["PySource1"]["Integral"]
  assert numpy.all(val == [2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0, 1002.0])

def test_picard_field1_oldint():
  val = stat["PicardProjection"]["Field1"]["OldIntegral"]
  assert numpy.all(abs(val - [1.33333333, 1.33333333, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0]) < 1.e-6)

def test_picard_source1_oldint():
  val = stat["PicardProjection"]["Source1"]["OldIntegral"]
  assert numpy.all(val == [2.0, 2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0])

def test_picard_pysource1_oldint():
  val = stat["PicardProjection"]["PySource1"]["OldIntegral"]
  assert numpy.all(val == [2.0, 2.0, 102.0, 202.0, 302.0, 402.0, 502.0, 602.0, 702.0, 802.0, 902.0])

def test_picard_field2_min_0():
  val = stat["PicardProjection"]["Field2"]["max"][0]
  assert numpy.all(val == [2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_field2_max_0():
  val = stat["PicardProjection"]["Field2"]["min"][0]
  assert numpy.all(val == [0.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_field2_int_0():
  val = stat["PicardProjection"]["Field2"]["Integral0"]
  assert numpy.all(abs(val - [2./3., 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0]) < 1.e-6)

def test_picard_source2_min_0():
  val = stat["PicardProjection"]["Source2"]["min"][0]
  assert numpy.all(val == [2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_source2_min_1():
  val = stat["PicardProjection"]["Source2"]["min"][1]
  assert numpy.all(val == [2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_source2_max_0():
  val = stat["PicardProjection"]["Source2"]["max"][0]
  assert numpy.all(val == [2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_source2_max_1():
  val = stat["PicardProjection"]["Source2"]["max"][1]
  assert numpy.all(val == [2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_source2_int_0():
  val = stat["PicardProjection"]["Source2"]["Integral0"]
  assert numpy.all(val == [2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_pysource2_int_0():
  val = stat["PicardProjection"]["PySource2"]["Integral0"]
  assert numpy.all(val == [2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0, 2002.0])

def test_picard_field2_oldint_0():
  val = stat["PicardProjection"]["Field2"]["OldIntegral0"]
  assert numpy.all(abs(val - [2./3., 2./3., 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0]) < 1.e-6)

def test_picard_source2_oldint_0():
  val = stat["PicardProjection"]["Source2"]["OldIntegral0"]
  assert numpy.all(val == [2.0, 2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0])

def test_picard_pysource2_oldint_0():
  val = stat["PicardProjection"]["PySource2"]["OldIntegral0"]
  assert numpy.all(val == [2.0, 2.0, 202.0, 402.0, 602.0, 802.0, 1002.0, 1202.0, 1402.0, 1602.0, 1802.0])

def test_picard_field2_min_1():
  val = stat["PicardProjection"]["Field2"]["max"][1]
  assert numpy.all(val == [2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_field2_max_1():
  val = stat["PicardProjection"]["Field2"]["min"][1]
  assert numpy.all(val == [0.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_field2_int_1():
  val = stat["PicardProjection"]["Field2"]["Integral1"]
  assert numpy.all(abs(val - [2./3., 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0]) < 1.e-6)

def test_picard_source2_int_1():
  val = stat["PicardProjection"]["Source2"]["Integral1"]
  assert numpy.all(val == [2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_pysource2_int_1():
  val = stat["PicardProjection"]["PySource2"]["Integral1"]
  assert numpy.all(val == [2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0, 3002.0])

def test_picard_field2_oldint_1():
  val = stat["PicardProjection"]["Field2"]["OldIntegral1"]
  assert numpy.all(abs(val - [2./3., 2./3., 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0]) < 1.e-6)

def test_picard_source2_oldint_1():
  val = stat["PicardProjection"]["Source2"]["OldIntegral1"]
  assert numpy.all(val == [2.0, 2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0])

def test_picard_pysource2_oldint_1():
  val = stat["PicardProjection"]["PySource2"]["OldIntegral1"]
  assert numpy.all(val == [2.0, 2.0, 302.0, 602.0, 902.0, 1202.0, 1502.0, 1802.0, 2102.0, 2402.0, 2702.0])

def test_snes_field1_min():
  val = stat["SNESProjection"]["Field1"]["min"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_max():
  val = stat["SNESProjection"]["Field1"]["max"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_int():
  val = stat["SNESProjection"]["Field1"]["Integral"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_source1_min():
  val = stat["SNESProjection"]["Source1"]["min"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_source1_max():
  val = stat["SNESProjection"]["Source1"]["max"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_source1_int():
  val = stat["SNESProjection"]["Source1"]["Integral"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_pysource1_int():
  val = stat["SNESProjection"]["PySource1"]["Integral"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_snes_field1_oldint():
  val = stat["SNESProjection"]["Field1"]["OldIntegral"]
  assert numpy.all(val == [0.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_snes_source1_oldint():
  val = stat["SNESProjection"]["Source1"]["OldIntegral"]
  assert numpy.all(val == [0.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_snes_pysource1_oldint():
  val = stat["SNESProjection"]["PySource1"]["OldIntegral"]
  assert numpy.all(val == [0.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_snes_field2_min_0():
  val = stat["SNESProjection"]["Field2"]["max"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_max_0():
  val = stat["SNESProjection"]["Field2"]["min"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_int_0():
  val = stat["SNESProjection"]["Field2"]["Integral0"]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_source2_min_0():
  val = stat["SNESProjection"]["Source2"]["min"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_source2_min_1():
  val = stat["SNESProjection"]["Source2"]["min"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_source2_max_0():
  val = stat["SNESProjection"]["Source2"]["max"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_source2_max_1():
  val = stat["SNESProjection"]["Source2"]["max"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_source2_int_0():
  val = stat["SNESProjection"]["Source2"]["Integral0"]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_pysource2_int_0():
  val = stat["SNESProjection"]["PySource2"]["Integral0"]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_snes_field2_oldint_0():
  val = stat["SNESProjection"]["Field2"]["OldIntegral0"]
  assert numpy.all(val == [0.0, 0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_snes_source2_oldint_0():
  val = stat["SNESProjection"]["Source2"]["OldIntegral0"]
  assert numpy.all(val == [0.0, 0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_snes_pysource2_oldint_0():
  val = stat["SNESProjection"]["PySource2"]["OldIntegral0"]
  assert numpy.all(val == [0.0, 0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_snes_field2_min_1():
  val = stat["SNESProjection"]["Field2"]["max"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_max_1():
  val = stat["SNESProjection"]["Field2"]["min"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_int_1():
  val = stat["SNESProjection"]["Field2"]["Integral1"]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_source2_int_1():
  val = stat["SNESProjection"]["Source2"]["Integral1"]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_pysource2_int_1():
  val = stat["SNESProjection"]["PySource2"]["Integral1"]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_snes_field2_oldint_1():
  val = stat["SNESProjection"]["Field2"]["OldIntegral1"]
  assert numpy.all(val == [0.0, 0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_snes_source2_oldint_1():
  val = stat["SNESProjection"]["Source2"]["OldIntegral1"]
  assert numpy.all(val == [0.0, 0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_snes_pysource2_oldint_1():
  val = stat["SNESProjection"]["PySource2"]["OldIntegral1"]
  assert numpy.all(val == [0.0, 0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

def test_picard_fieldeval1_min():
  val = stat["PicardProjection"]["FieldEval1"]["min"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_fieldeval1_max():
  val = stat["PicardProjection"]["FieldEval1"]["max"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_fieldeval1_int():
  val = stat["PicardProjection"]["FieldEval1"]["Integral"]
  assert numpy.all(val == [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

def test_picard_fieldeval1_oldint():
  val = stat["PicardProjection"]["FieldEval1"]["OldIntegral"]
  assert numpy.all(val == [0.0, 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0])

def test_picard_fieldeval2_min_0():
  val = stat["PicardProjection"]["FieldEval2"]["max"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_fieldeval2_max_0():
  val = stat["PicardProjection"]["FieldEval2"]["min"][0]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_fieldeval2_int_0():
  val = stat["PicardProjection"]["FieldEval2"]["Integral0"]
  assert numpy.all(val == [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0])

def test_picard_fieldeval2_oldint_0():
  val = stat["PicardProjection"]["FieldEval2"]["OldIntegral0"]
  assert numpy.all(val == [0.0, 0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])

def test_picard_fieldeval2_min_1():
  val = stat["PicardProjection"]["FieldEval2"]["max"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_fieldeval2_max_1():
  val = stat["PicardProjection"]["FieldEval2"]["min"][1]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_fieldeval2_int_1():
  val = stat["PicardProjection"]["FieldEval2"]["Integral1"]
  assert numpy.all(val == [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0])

def test_picard_fieldeval2_oldint_1():
  val = stat["PicardProjection"]["FieldEval2"]["OldIntegral1"]
  assert numpy.all(val == [0.0, 0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0])

