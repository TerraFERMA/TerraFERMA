
from buckettools.statfile import parser
import numpy

stat = parser("adaptdt.stat")

def test_mincourant():
  val = stat["CourantNumber"]["AltCourantNumber"]["min"][1:].min()
  assert numpy.all(val==10.)

def test_maxcourant():
  val = stat["CourantNumber"]["AltCourantNumber"]["max"][1:].max()
  assert numpy.all(val==15.0)

def test_dt():
  val = stat["dt"]["value"]
  expectedval = numpy.array([ 1.    ,  1.    ,  1.    , \
                              1.1   ,  1.1   ,  \
                              1.21  ,  1.21  ,  \
                              1.331 ,  1.331 ,  \
                              1.4641,  1.4641,  \
                              1.5   ,  1.5   ,  1.5   ,    1.5   ,  1.5   ,  1.5   ])
  assert numpy.all(abs(val-expectedval) < 1.e-14)

