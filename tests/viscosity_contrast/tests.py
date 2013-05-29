
from buckettools.statfile import parser
import numpy

conv = parser("dc_Benton_Solver_ksp.conv")
logfile = file("bucket.log-0", 'r')
loglines = logfile.readlines()
logfile.close()

convit = conv["KSPIteration"]['value'].max()
logits = numpy.array([int(line.split()[0]) for line in loglines if line.find("KSP Residual norm") > 0 ])
logmaxits = logits[logits[1:]-logits[:-1] < 0]

def test_conv_it():
  assert abs(convit-36) <=2

def test_field_max():
  assert numpy.all(logmaxits==convit)

