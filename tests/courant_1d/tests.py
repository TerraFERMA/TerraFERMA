
from buckettools.statfile import parser
import numpy

stat = parser("courant.stat")

def test_mincourant():
  val = stat["CourantNumber"]["CourantNumber"]["min"][1:]
  assert numpy.all(val==10.)

def test_maxcourant():
  val = stat["CourantNumber"]["CourantNumber"]["max"][1:]
  assert numpy.all(val==10.)

