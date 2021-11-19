"""
File: magmasinc.py

Descirption: Compute the d-dimensional solitary waves for several magma equations
solved using sinc nonlinear collocation.

Solves all cases, with special routines for:

n=3, m=0:  solwave_mck
n=2, m=1:  solwave_con
n>1, m=1:  solwave_m1
n>1, m<>1: solwave_gen

Author: Gideon Simpson, simpson@math.toronto.edu
"""

from . import sinc_eo, magma1d
import numpy as np
from numpy import diag, dot, sqrt, pi, linspace, log, cosh, max
from scipy.optimize.minpack import fsolve
from scipy.linalg import solve

def solwave_mck_eq(u, params):
	"""
	The solitary wave PDE posed as a nonlinear collocation problem
	in the case of n=3, m=0. 
	
	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices
	
	Outputs: a vector of the residual values associated with the
	discretized equation
	"""
	
	# Extract equation parameters
	c = params['c']
	d = params['d']

	# Extract matrices 
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	return -u + (1. / c) * (f**3 -	  1.) + (f**3) * dot(D2, u) \
	       + (d - 1.) * dot(IN, (f**3) * dot(D1, dot(D1_r, u)))

def solwave_mck_eq_jac(u, params):
	"""
	The Jacobian of the solitary wave PDE posed as a nonlinear
	collocation problem in the case of n=3, m=0.

	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices
	
	Outputs: the Jacobian matrix
	"""

	# Extract equation parameters
	c = params['c']
	d = params['d']

	# Extract matrices 
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	
	return diag(-1. + (3. / c) * f**2 + dot(D2, u) * 3. * f**2) \
	       + dot(diag(f**3), D2) \
	       + (d - 1.) * dot(IN, dot(diag(f**3), dot(D1, D1_r))) \
	       + (d - 1.) * dot(IN, dot(diag(3.* f**2), diag(dot(D1, dot(D1_r, u)))))

def solwave_mck(c, d, M, dist = pi / 2., xtol = 1.e-12, verbose = False):
	"""
	Compute the collocation points and values for the d-dimensional
	McKenzie solitary wave equation.  n=3, m=0.

	Inputs: c - soliton parameter, 

	d - dimension, 

	M - No. of points, 

	[dist = pi/2] - distance from the x axis to the nearest
	singularity in the complex plane, 

	[xtol = 1.e-12] - solver tolerance,

	[verbose=True/False]

	Outputs: r, f - collocation points and the soliton
	"""
		
	decay = sqrt(1. - 3. / c)

	h = sqrt((pi * dist) / (decay * M))
	
	r = sinc_eo.points(M, h)
	
	u0 = magma1d.solwave_mck(r, c) - 1.

	d0 = 1.
	
	D2 = sinc_eo.D2_e(M, h)
	D1 = sinc_eo.D1_eo(M, h)
	D1_r = sinc_eo.D1_x_e(M, h)
	IN1 = sinc_eo.IN1_oe(M, h)
	
	params = {'c': c, 'd': d0, 
		  'D2': D2, 'D1':D1, 'D1_r':D1_r, 'IN': IN1}

	if d < 2.:
		
		params['d'] = d

		u = fsolve(lambda v: solwave_mck_eq(v, params), u0, 
			   fprime = lambda v:solwave_mck_eq_jac(v, params), 
			   xtol= xtol)
	
	else:
		# Initial value for the delta dimension value.  
		# This can be tuned as needed
		delta_dim = (d - 1.)/10.

		# Loop will terminate if the delta dim value gets too
		# small, this may be tunded as needed
		delta_dim_min = 1.e-6

		d0 = 1.

		while d0 < d and delta_dim > delta_dim_min:

			d1 = min([d0 + delta_dim, d])
			params['d'] = d1

			if verbose:
				print(' computing with d = ', d1)
		
			u1, infodict, ier, mesg = \
			    fsolve(lambda v:solwave_mck_eq(v,params), u0, \
				   fprime =lambda v:solwave_mck_eq_jac(v, params), \
				   xtol= xtol,
				   full_output=1)



			if ier == 1 and np.max(np.abs(u1)) > 1.e-8:

				d0 = d1
				u0 = u1

			else:
				if verbose:
					if np.max(np.abs(u1)) <= 1.e-8:
						print(" converged to the zero solution, adjusteing delta_dim")
					else:
						print(" solver failed to converge, adjusting delta_dim")
						print(" solver err = ", ier)

				delta_dim = delta_dim / 2.
			if verbose:
				jac = solwave_mck_eq_jac(u1, params)
				condnum = np.linalg.cond(jac)
				print(" condition number = ", condnum)

				
		u = u0

	if d0 < d:
		print(" Failed to converge to the soliton solution, last value of d =", d0)
		f = 0.0
	else:
		f = 1. + u

	return r, f

def solwave_con_eq(u, params):
	"""
	The solitary wave PDE posed as a nonlinear collocation problem
	in the case of n=2, m=1.

	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices

	Outputs: a vector of the residual values associated with the
	discretized equation	
	"""

	# Extract equation parameters	
	c = params['c']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u


	return -u + (1. / c) * (f**2 - 1.) + f**2 * dot(D2, log(f) ) \
	       + (d - 1.) * dot(IN, f**2 * dot(D1, dot(D1_r,log(f) ) ) )

def solwave_con_eq_jac(u, params):
	"""
	The Jacobian of the solitary wave PDE posed as a nonlinear
	collocation problem in the case of n=2, m=1.

	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices

	Outputs: the Jacobian matrix
	"""

	# Extract equation parameters	
	c = params['c']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u

	return  diag(-1. + (2. / c) * f +  2. * f * dot(D2, log(f))  ) \
	       + dot(diag(f**2), dot(D2, diag(1./f)) ) \
	       + (d - 1.) *  dot(IN, dot(diag(f**2), dot(D1, dot(D1_r, diag(1./f))))) \
	       + (d - 1.) *  dot(IN, dot(diag(2. * f), diag( dot(D1, dot(D1_r, log(f))) ) ) ) 


def solwave_con(c, d, M, dist = pi / 2., xtol = 1.e-12, verbose = False):
	"""
	Compute the collocation points and values for the d-dimensional
	conduit solitary wave equation.  n=2, m=1.

	Inputs: c - soliton parameter, 

	d - dimension, 

	M - No. of points, 

	[dist = pi/2] - distance from the x axis to the nearest
	singularity in the complex plane, 

	[xtol = 1.e-12] - solver tolerance,

	[verbose=True/False]

	Outputs: r, f - collocation points and the soliton
	"""

	# Iterate in soliton parameter c, if neccessary, to achieve a good
	# guess.  This may be replaced with another algorithm if a better
	# initial guess for the 1D soliton is available.
	if verbose:
		print(" iterating in c")

	# Initial value for delta c.  This can be tuned as needed.
	delta_c = .5 * 2./c

	# Loop will terminate if the delta c value gets too small, this may be
	# adjusted as needed
	delta_c_min = 1.e-6

	c0 = 2.
	d0 = 1.
	
	while c0 < c and delta_c > delta_c_min:

		c1 = min([c0 + delta_c, c])

		decay = sqrt(1. - 2. / c1)

		h = sqrt((pi * dist) / (decay * M))

		r1 = sinc_eo.points(M, h)

		D2 = sinc_eo.D2_e(M, h)
		D1 = sinc_eo.D1_eo(M, h)
		D1_r = sinc_eo.D1_x_e(M, h)
		IN1 = sinc_eo.IN1_oe(M, h)

		params = {'c': c1, 'd': d0,
			  'D2': D2, 'D1':D1, 'D1_r':D1_r, 'IN': IN1}

		try:
			uguess = sinc_eo.sincinterp_e(r0,u0,r1)
		except:
			uguess = 3. * (1. - 2. / c1) / cosh(.5 * decay * r1)**2

		if verbose:
			print(' computing with c = ', c1)

		u1, infodict, ier, mesg = \
		    fsolve(lambda v:solwave_con_eq(v, params), uguess, \
			   fprime = lambda v:solwave_con_eq_jac(v,params), \
			   xtol= xtol,
			   full_output=1)

		if ier == 1 and np.max(np.abs(u1)) > 1.e-10:

			c0 = c1
			u0 = u1
			r0 = r1

		else:
			if verbose:
				if np.max(np.abs(u1)) <= 1.e-8:
					print(" converged to the zero solution, adjusteing delta_c")
				else:
					print(" solver failed to converge, adjusting delta_c")
					print(" solver err = ", ier)

			delta_c = delta_c / 2.
			
		if verbose:
			jac = solwave_con_eq_jac(u1, params)
			condnum = np.linalg.cond(jac)
			print(" condition number = ", condnum)
			

	u = u0
	r = r0

	# Iterate in dimension, if neccessary	
	if d > 1.:
		
		if verbose:
			print(" iterating in d")

		delta_dim = (d - 1.) / 10.

		# Loop will terminate if the delta dim value gets too small, this may be
		# adjusted as needed
		delta_dim_min = 1.e-6

		d0 = 1.

		u0 = u

		while d0 < d and delta_dim > delta_dim_min:

			d1 = min([d0 + delta_dim, d])
			params['d'] = d1

			if verbose:
				print(' computing with d = ', d1)
		
			u1, infodict, ier, mesg = \
			    fsolve(lambda v:solwave_con_eq(v,params), u0, \
				   fprime = lambda v:solwave_con_eq_jac(v,params), \
				   xtol= xtol,
				   full_output=1)

			if ier == 1 and np.max(np.abs(u1)) > 1.e-8:

				d0 = d1
				u0 = u1

			else:
				if verbose:
					if np.max(np.abs(u1)) <= 1.e-8:
						print(" converged to the zero solution, adjusteing delta_dim")
					else:
						print(" solver failed to converge, adjusting delta_dim")
						print(" solver err = ", ier)
						print(mesg)

				delta_dim = delta_dim / 2.

			if verbose:
				jac = solwave_con_eq_jac(u1, params)
				condnum = np.linalg.cond(jac)
				print(" condition number = ", condnum)
				
		u = u0

	if d0 < d:

		print(" Failed to converge to the soliton solution, last value of d =", d0)

		f = 0.0 * r
		return r, f

	else:
		f = 1. + u

	return r, f


def solwave_gen_eq(u, params):
	"""
	The solitary wave PDE posed as a nonlinear collocation problem
	in the general case of n>1, and m anything but 1.
	
	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices
	
	Outputs: a vector of the residual values associated with the
	discretized equation
	"""
	
	# Extract equation parameters	
	c = params['c']
	n = params['n']
	m = params['m']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	return -u + (1. / c) * (f**n - 1.) \
	       + (1. / (1. - m)) * (f**n) * dot(D2, (f**(1. - m) - 1.)) \
	       + (d - 1.) / (1. - m) \
	       * dot(IN, (f**n) * dot(D1, dot(D1_r, (f**(1. - m) - 1.))))


def solwave_gen_eq_jac(u, params):
	"""
	The Jacobian of the solitary wave PDE posed as a nonlinear
	collocation problem in the case of n>1, m m anything but 1.

	Inputs: u - the soliton, 
	params - data structure containing the equation parameters and
	the collocation matrices
	
	Outputs: the Jacobian matrix
	"""

	# Extract equation parameters	
	c = params['c']
	n = params['n']
	m = params['m']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	
	return diag(-1. + (n / c) * f**(n - 1.) \
		      + (1. / (1. - m)) * n * f**(n-1.) * dot(D2,(f**(1. - m) - 1.))) \
		+ dot( diag(f**n), dot(D2,diag(f**(-m)))) \
		+ (d - 1.) * dot(IN, dot(diag(f**n), dot(D1,dot(D1_r, diag(f**(-m)) )))) \
		+ (d - 1.) / (1. - m) \
		* dot(IN, dot(diag(n * f**(n - 1.)), diag(dot(D1, dot(D1_r, (f**(1.-m) - 1.))))))

def solwave_gen(c, n, m, d, M, dist = pi / 2., xtol = 1.e-12, verbose = False):
	"""
	Compute the collocation points and values for the general d-dimensional
	solitary wave equation with n > 1, and m not equal to 1

	Inputs: c - soliton parameter, 

	n - first nonlinearity, 

	m - second nonlinearity,

	d - dimension, 

	M - No. of points, 

	[dist = pi/2] - distance from the x axis to the nearest
	singularity in the complex plane, 

	[xtol = 1.e-12] - solver tolerance,

	[verbose=True/False]


	Outputs: r, f - collocation points and the soliton
	"""

	# Iterate in soliton parameter, if neccessary, to achieve a good
	# guess.  This may be replaced with another algorithm if a better
	# initial guess for the 1D soliton is available.
	if verbose:
		print(" iterating in c")

	# Initial value for delta c.  This can be tuned as needed.
	delta_c = .5 * n / c

	# Loop will terminate if the delta c value gets too small, this may be
	# adjusted as needed
	delta_c_min = 1.e-6

	c0 = n
	d0 = 1.

	while c0 < c and delta_c > delta_c_min:

		c1 = min([c0 + delta_c, c])

		decay = sqrt(1. - n / c1)

		h = sqrt((pi * dist) / (decay * M))

		r1 = sinc_eo.points(M, h)

		D2 = sinc_eo.D2_e(M, h)
		D1 = sinc_eo.D1_eo(M, h)
		D1_r = sinc_eo.D1_x_e(M, h)
		IN1 = sinc_eo.IN1_oe(M, h)

		params = {'c': c1, 'n': n, 'm': m, 'd': d0,
			  'D2': D2, 'D1': D1, 'D1_r': D1_r, 'IN': IN1}

		try:
			uguess = sinc_eo.sincinterp_e(r0, u0, r1)
		except:
			uguess = 3. * (1. - n / c1) / cosh(.5 * decay * r1)**2

		if verbose:
			print(' computing with c = ', c1)

		u1, infodict, ier, mesg = \
		    fsolve(lambda v: solwave_gen_eq(v,params), uguess, \
			   fprime = lambda v:solwave_gen_eq_jac(v, params), 
			   xtol= xtol, 
			   full_output=1)

		if ier == 1 and np.max(np.abs(u1)) > 1.e-10:

			c0 = c1
			u0 = u1
			r0 = r1

		else:
			if verbose:
				if np.max(np.abs(u1)) <= 1.e-8:
					print(" converged to the zero solution, adjusteing delta_c")
				else:
					print(" solver failed to converge, adjusting delta_c")
					print(" solver err = ", ier)

			delta_c = delta_c / 2.

	u = u0
	r = r0

	# Iterate in dimension, if neccessary
	if d > 1.:
		if verbose:
			print(" iterating in d")

		delta_dim = (n - 1.) / 10.

		# Loop will terminate if the delta dim value gets too
		# small, this can be adjusted as desired 
		delta_dim_min = 1.e-6

		d0 = 1.

		u0 = u

		while d0 < d and delta_dim > delta_dim_min:

			d1 = min([d0 + delta_dim, d])

			params['d'] = d1

			if verbose:
				print(' computing with d = ', d1)
		
			u1, infodict, ier, mesg = \
			    fsolve(lambda v: solwave_gen_eq(v, params), u0, \
				   fprime = lambda v:solwave_gen_eq_jac(v, params), 
				   xtol= xtol,
				   full_output=1)


			if ier == 1 and np.max(np.abs(u1)) > 1.e-8:

				d0 = d1
				u0 = u1

			else:
				if verbose:
					if np.max(np.abs(u1)) <= 1.e-8:
						print(" converged to the zero solution, adjusteing delta_dim")
					else:
						print(" solver failed to converge, adjusting delta_dim")
						print(" solver err = ", ier)

				delta_dim = delta_dim / 2.
				
		u = u0
		
	if d0 < d:
		print(" Failed to converge to the soliton solution, last value of d =", d0)
		f = 0.0
	else:
		f = 1. + u

	return r, f

def solwave_m1_eq(u, params):
	"""
	The solitary wave PDE posed as a nonlinear collocation problem
	in the general case of n>1 and m = 1.
	
	Inputs: u - the soliton, c - the soliton parameter, d - the
	dimension, mats - the various collocation matrices
	
	Outputs: a vector of the residual values associated with the
	discretized equation
	"""
	
	# Extract equation parameters	
	c = params['c']
	n = params['n']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	return -u + (1. / c) * (f**n - 1.) + (f**n) * dot(D2, log(f)) \
	       + (d - 1.) * dot(IN, (f**n) * dot(D1, dot(D1_r, log(f) ) ) ) 

def solwave_m1_eq_jac(u, params):
	"""
	The Jacobian of the solitary wave PDE posed as a nonlinear
	collocation problem in the case of n>1 and m = 1.
	
	Inputs: u - the soliton, c - the soliton parameter, d - the
	dimension, mats - the various collocation matrices
	
	Outputs: a vector of the residual values associated with the
	discretized equation
	"""
	
	# Extract equation parameters	
	c = params['c']
	n = params['n']
	d = params['d']

	# Extract matrices		
	D1 = params['D1']
	D1_r = params['D1_r']
	D2 = params['D2']
	IN = params['IN']

	f = 1. + u
	return diag(-1. + n / c * f**(n - 1.) + n * f**(n - 1.) * dot(D2, log(f))) \
	       + dot(diag(f**n), dot(D2, diag(1./f)) ) \
	       + (d - 1.) *  dot(IN, dot(diag(f**n), dot(D1, dot(D1_r, diag(1./f))))) \
	       + (d - 1.) *  dot(IN, dot(diag(n * f**(n-1.)), diag(dot(D1, dot(D1_r, log(f))  ) )))

def solwave_m1(c, n, d, M, dist = pi / 2., xtol = 1.e-12, verbose = False):
	"""
	Compute the collocation points and values for the general d-dimensional
	solitary wave equation with n > 1, and m=1

	Inputs: c - soliton parameter, 

	n - first nonlinearity, 

	d - dimension, 

	M - No. of points, 

	[dist = pi/2] - distance from the x axis to the nearest
	singularity in the complex plane, 

	[xtol = 1.e-12] - solver tolerance,

	[verbose=True/False]

	Outputs: r, f - collocation points and the soliton
	"""

	# Iterate in soliton parameter, if neccessary, to achieve a good
	# guess.  This may be replaced with another algorithm if a better
	# initial guess for the 1D soliton is available.
	if verbose:
		print(" iterating in c")

	# Initial value for delta c.  This can be tuned as needed.
	delta_c = .5 * n / c

	# Loop will terminate if the delta c value gets too small, this may be
	# adjusted as needed
	delta_c_min = 1.e-6

	c0 = n
	d0 = 1.

	while c0 < c and delta_c > delta_c_min:

		c1 = min([c0 + delta_c, c])

		decay = sqrt(1. - n / c1)

		h = sqrt((pi * dist) / (decay * M))

		r1 = sinc_eo.points(M, h)

		D2 = sinc_eo.D2_e(M, h)
		D1 = sinc_eo.D1_eo(M, h)
		D1_r = sinc_eo.D1_x_e(M, h)
		IN1 = sinc_eo.IN1_oe(M, h)

		params = {'c': c1, 'n': n, 'd': d0, 
			  'D2': D2, 'D1': D1, 'D1_r': D1_r, 'IN': IN1}

		try:
			uguess = sinc_eo.sincinterp_e(r0,u0,r1)
		except:
			uguess = 3. * (1. - n / c1) / cosh(.5 * decay * r1)**2

		if verbose:
			print(' computing with c = ', c1)

		u1, infodict, ier, mesg = \
		    fsolve(lambda v: solwave_m1_eq(v,params), uguess, \
#			   fprime = lambda v:solwave_m1_eq_jac(v, params), \
			   xtol= xtol,
			   full_output=1)

		if ier == 1 and np.max(np.abs(u1)) > 1.e-10:

			c0 = c1
			u0 = u1
			r0 = r1

		else:
			if verbose:
				if np.max(np.abs(u1)) <= 1.e-8:
					print(" converged to the zero solution, adjusteing delta_c")
				else:
					print(" solver failed to converge, adjusting delta_c")
					print(" solver err = ", ier)

			delta_c = delta_c / 2.

	u = u0
	r = r0

	# Iterate in dimension, if neccessary
	if d > 1.:

		if verbose:
			print(" iterating in d")

		delta_dim = (n - 1.) / 10.

		# Loop will terminate if the delta dim value gets too
		# small, this can be adjusted as desired 
		delta_dim_min = 1.e-6

		d0 = 1.

		u0 = u

		while d0 < d and delta_dim > delta_dim_min:

			d1 = min([d0 + delta_dim, d])

			params['d'] = d1

			if verbose:
				print(' computing with d = ', d1)
		
			u1, infodict, ier, mesg = \
			    fsolve(lambda v: solwave_m1_eq(v, params), u0, \
#				   fprime = lambda v:solwave_m1_eq_jac(v, params),
				   xtol= xtol,
				   full_output=1)


			if ier == 1 and np.max(np.abs(u1)) > 1.e-8:

				d0 = d1
				u0 = u1

			else:
				if verbose:
					if np.max(np.abs(u1)) <= 1.e-8:
						print(" converged to the zero solution, adjusteing delta_dim")
					else:
						print(" solver failed to converge, adjusting delta_dim")
						print(" solver err = ", ier)

				delta_dim = delta_dim / 2.
				
		u = u0
		
	if d0 < d:
		print(" Failed to converge to the soliton solution, last value of d =", d0)
		f = 0.0
	else:
		f = 1. + u

	return r, f

