
from buckettools.statfile import parser
import numpy

stat = parser("poisson.stat")

def test_integral():
  val = stat["SNESPoisson"]["Field"]["Integral"][1]
  assert abs(val)<1.e-10

