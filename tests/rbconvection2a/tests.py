
from buckettools.statfile import parser
import numpy
from math import sqrt

stat = parser("rbconvection.stat")

#def test_timestepcount():
#  val = stat["timestep"]["value"][-1]
#  assert abs(val - 2167) < 50
#
#def test_elapsedtime():
#  val = stat["ElapsedTime"]["value"][-1]
#  assert abs(val - 2167.0) < 5.e1
#
def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])
  assert abs(val - 471.1922e-4) < 0.005

def test_nu():
  val = -1.0*(stat["Stokes"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 10.1565 < 0.5

