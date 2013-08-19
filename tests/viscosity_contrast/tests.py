# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


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

