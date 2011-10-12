
from buckettools.statfile import parser
import numpy
from math import sqrt

det = parser("subduction.det")

def test_T_11_11():
  val = det["Temperature"]["Temperature"]["SlabPoint"][0,-1]-273.
  assert abs(val - 587.0) < 1.0

def test_T_slab():
  val = sqrt(sum((det["Temperature"]["Temperature"]["Slab"][:,-1]-273.)**2)/36.)
  assert abs(val - 609.0) < 1.0

def test_T_wedge():
  val = sqrt(sum((det["Temperature"]["Temperature"]["Wedge"][:,-1]-273.)**2)/78.)
  assert abs(val - 1005.0) < 1.0

