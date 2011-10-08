
from buckettools.statfile import parser
import numpy

stat = parser("region_ids.stat")
det = parser("region_ids.det")

def test_field_min():
  val = stat["Regions"]["Field"]["min"]
  assert numpy.all(val==32.0)

def test_field_max():
  val = stat["Regions"]["Field"]["max"]
  assert numpy.all(val==35.0)

def test_coeff_min():
  val = stat["Regions"]["Coefficient"]["min"]
  assert numpy.all(val==32.0)

def test_coeff_max():
  val = stat["Regions"]["Coefficient"]["max"]
  assert numpy.all(val==35.0)

def test_field_lower():
  val = det["Regions"]["Field"]["Lower"]
  assert numpy.all(val==35.0)

def test_field_lowermiddle():
  val = det["Regions"]["Field"]["LowerMiddle"]
  assert numpy.all(val==34.0)

def test_field_uppermiddle():
  val = det["Regions"]["Field"]["UpperMiddle"]
  assert numpy.all(val==33.0)

def test_field_upper():
  val = det["Regions"]["Field"]["Upper"]
  assert numpy.all(val==32.0)

def test_coeff_lower():
  val = det["Regions"]["Coefficient"]["Lower"]
  assert numpy.all(val==35.0)

def test_coeff_lowermiddle():
  val = det["Regions"]["Coefficient"]["LowerMiddle"]
  assert numpy.all(val==34.0)

def test_coeff_uppermiddle():
  val = det["Regions"]["Coefficient"]["UpperMiddle"]
  assert numpy.all(val==33.0)

def test_coeff_upper():
  val = det["Regions"]["Coefficient"]["Upper"]
  assert numpy.all(val==32.0)

