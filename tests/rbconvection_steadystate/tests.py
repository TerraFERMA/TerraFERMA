
from buckettools.statfile import parser
import numpy
from math import sqrt

stat = parser("rbconvection.stat")

def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2Norm"][-1])
  assert abs(val - 42.865) < 0.01

