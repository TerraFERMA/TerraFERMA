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

"""
One dimensional solitary waves of the n=3, m=0 magma equation
of D. McKenzie.  Makes use of the first integral from
Barcilon & Richter.

Gideon Simpson
simpson@math.toronto.edu
"""

import numpy
import scipy
from numpy import sqrt, log
from scipy.optimize.zeros import brentq


def solwave_amp_mck(c):
	"""
	Compute solitary wave amplitude for speed c.
	"""
	
	return .5 * (c-1.)

def solwave_implicit_mck(f, x, amp):
	"""
	Evaluate implicit solitary wave equation.
	"""
	
	return x - sqrt(amp + .5) * (2. * sqrt(amp - f) - log((sqrt(amp - 1.) - 
		sqrt(amp - f))/(sqrt(amp - 1.) + sqrt(amp - f)))/sqrt(amp - 1.))


def solwave_mck(x, c):
	"""
	Return solitary wave values at x, f_c(x).
	"""
	
	amp = solwave_amp_mck(c)

	def bracket(xx):
		"""
		Construct bracket for root solving about xx.
		"""
		
		aa = .5 * (1. + amp)
		bb = amp
		
		while (solwave_implicit_mck(aa, xx, amp)>0 and 
			numpy.isfinite(solwave_implicit_mck(aa, xx, amp))):
			bb = aa
			aa = .5 * (aa + 1.)
		
		if numpy.isfinite(solwave_implicit_mck(aa, xx, amp)):
			return aa, bb
		else:
			return -1., -1.

	
	if numpy.isscalar(x):
		
		a, b = bracket(x)
		if a > 0. and b > 0.:
			f = brentq(solwave_implicit_mck, a, b, args = (x, amp))
		else:
			f = 1.
	
	else:	
		
		f = numpy.zeros(len(x))
		
		for j in range(0, len(x)):
			
			a, b = bracket(x[j])
			if a > 0. and b > 0.:
				f[j] = brentq(solwave_implicit_mck, a, b, args = (x[j], amp))
			else:
				f[j] = 1.
	
	return f
	
	