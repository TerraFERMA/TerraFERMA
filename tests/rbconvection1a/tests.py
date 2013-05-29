
from buckettools.statfile import parser
import numpy
from math import sqrt

stat = parser("rbconvection.stat")
det = parser("rbconvection.det")
steady = parser("rbconvection.steady")

def test_timestepcount():
  val = stat["timestep"]["value"][-1]
  assert abs(val - 30) < 1

def test_elapsedtime():
  val = stat["ElapsedTime"]["value"][-1]
  assert abs(val - 1640.9259654) < 5.e1

def test_v_rms():
  val = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])
  assert abs(val - 42.865e-4) < 0.01

def test_nu():
  val = -1.0*(stat["Stokes"]["Temperature"]["TopSurfaceIntegral"][-1])
  assert val - 4.9 < 0.05

def test_extremumloc():
  val = (det["Stokes"]["Temperature"]["Array"][0:129/2,-1]).argmin()*1./128.
  assert abs(val - 0.2265625) < 0.01

def test_extremum():
  val = (det["Stokes"]["Temperature"]["Array"][0:129/2,-1]).min()
  assert abs(val - 0.4222) < 0.01

def test_steady():
  val = max(steady["Stokes"]["Velocity"]["change(linf)"][-1], \
            steady["Stokes"]["Pressure"]["change(linf)"][-1], \
            steady["Stokes"]["Temperature"]["change(linf)"][-1])
  assert val < 1.e-5

def test_div():
  val = max(numpy.abs(stat["Divergence"]["Divergence"]["max"]).max(), \
            numpy.abs(stat["Divergence"]["Divergence"]["min"]).min())
  assert val < 1.e-6

