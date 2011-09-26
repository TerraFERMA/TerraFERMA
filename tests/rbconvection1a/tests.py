
from buckettools.statfile import parser
import numpy
from math import sqrt

stat = parser("rbconvection.stat")

def test_timestepcount():
  val = stat["timestep"]["value"][-1]
  assert abs(val - 66) < 1

def test_elapsedtime():
  val = stat["ElapsedTime"]["value"][-1]
  assert abs(val - 3.3e3) < 5.e1

def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])
  assert abs(val - 42.865e-4) < 0.01

def test_nu():
  val = -1.0*(stat["Stokes"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

