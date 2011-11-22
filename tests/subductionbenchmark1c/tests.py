
from buckettools.statfile import parser
import numpy
from math import sqrt

det = parser("subduction.det")

def test_T_11_11():
  val = det["Solid"]["Temperature"]["SlabPoint"][0,-1]-273.
  test = 388.5
  assert abs(val - test) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_slab():
  val = sqrt(sum((det["Solid"]["Temperature"]["Slab"][:,-1]-273.)**2)/36.)
  test = 504.0
  assert abs(val - test) < 1.0
  print "\tvalue=", val, "test=", test, " ",

def test_T_wedge():
  val = sqrt(sum((det["Solid"]["Temperature"]["Wedge"][:,-1]-273.)**2)/78.)
  test = 853.5
  assert abs(val - test) < 1.0
  print "\tvalue=", val, "test=", test, " ",

