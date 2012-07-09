
from buckettools.statfile import parser
import numpy
from math import sqrt

stat = parser("rbconvection.stat")

def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2Norm"][-1])
  assert abs(val - 42.865) < 0.01

def test_nu():
  val = -1.0*(stat["Temperature"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

def test_v2_rms():
  val = sqrt(stat["Stokes2"]["Velocity"]["L2Norm"][-1])
  assert abs(val - 42.865) < 0.01

def test_nu2():
  val = -1.0*(stat["Temperature2"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

def test_v3_rms():
  val = sqrt(stat["Stokes3"]["Velocity"]["L2Norm"][-1])
  assert abs(val - 42.865) < 0.01

def test_nu3():
  val = -1.0*(stat["Temperature3"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

