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
	
	