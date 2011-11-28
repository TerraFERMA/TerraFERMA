
from buckettools.statfile import parser
import numpy
from math import sqrt

det = parser("subduction.det")

def test_T_11_11():
  val = det["Solid"]["Temperature"]["SlabPoint"][0,-1]-273.
  test = 587.0
  assert abs(val - 587.0) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_slab():
  val = sqrt(sum((det["Solid"]["Temperature"]["Slab"][:,-1]-273.)**2)/36.)
  test = 609.0
  assert abs(val - 609.0) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_wedge():
  val = sqrt(sum((det["Solid"]["Temperature"]["Wedge"][:,-1]-273.)**2)/78.)
  test = 1000.5
  assert abs(val - 1005.0) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_runtime():
  val = det["ElapsedWallTime"]["value"][-1]
  test = 20.
  assert ( val <= test)
  print "\tElapsed Time= ",val," s ",
