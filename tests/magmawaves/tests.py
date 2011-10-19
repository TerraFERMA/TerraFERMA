from buckettools.statfile import parser
import numpy
from math import sqrt

det = parser("magmawaves.det")
stat = parser("magmawaves.stat")

def test_Porosity_center():
  val = det["magma"]["Porosity"]["x0"][0,-1]
  test = 2.3359664311e+00
  assert abs(val - test) < 0.001
  print "\tvalue=", val, "test=", test, " ",

def test_Pressure_center():
  val = det["magma"]["Pressure"]["x0"][0,-1]
  test = 0
  assert abs(val - test) < 1.e-5
  print "\tvalue=", val, "test=", test, " ",

def test_Porosity_integral():
  val = stat["magma"]["Porosity"]["Integral"][-1]
  test = 1.0378229285e+00
  assert abs(val - test) < 0.00001
  print "\tvalue=", val, "test=", test, " ",

def test_Pressure_integral():
  val = stat["magma"]["Pressure"]["IntegralPressure"][-1]
  test = 0.
  assert abs(val - test) < 1.e-10
  print "\tvalue=", val, "test=", test, " ",

def test_timestep():
  val = stat["timestep"]["value"][-1]
  test = 1282
  assert abs(val - test) < 1
  print "\tvalue=", val, "test=", test, " ",

def test_timestep():
  val = stat["timestep"]["value"][-1]
  test = 1282
  assert abs(val - test) < 1
  print "\tvalue=", val, "test=", test, " ",

def test_walltime():
  val = stat["ElapsedWallTime"]["value"][-1]
  test = 622.
  assert abs(val - test) < 30.
  print "\tvalue=", val, "test=", test, " ",
