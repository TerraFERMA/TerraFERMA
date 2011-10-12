
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
