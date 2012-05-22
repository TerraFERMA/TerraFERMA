
from buckettools.statfile import parser
import numpy

stat = parser("advection.stat")

def test_upper_bound():
  val = stat["System"]["Field"]["max"]
  assert numpy.all(val<=1.0)

def test_lower_bound():
  val = stat["System"]["Field"]["min"]
  assert numpy.all(val>=0.0)

