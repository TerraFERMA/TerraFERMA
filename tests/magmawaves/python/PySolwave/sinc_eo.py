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
Matrices and collocation points for using sinc spectral methods
on problems with even and odd symmetry.

Gideon Simpson
simpson@math.toronto.edu
"""

import sincmats
from numpy import dot, sqrt, vstack, diag, hstack, zeros, fliplr, arange, pi, sinc
#from scipy.interpolate.interpolate import spline


def points(N, h):
	"""
	Compute the collocation points 0, h, 2h, ..., Nh
	"""
	
	return sincmats.points(2 * N + 1, h)[-N-1:]

def dot_eo(u, v):
	"""
	Dot product of two even or two odd vectors.
	"""
	
	return 2*dot(u,v) - u[0]*v[0]

def norm_eo(u, h):
	"""
	l2 norm of an even or odd vector, with h weighting.
	"""
	
	return sqrt(h * evendot(u, u))

def D2_e(N, h):
	"""
	Compute the d2 matrix for even functions.
	"""
	
	DM = sincmats.deriv_mats(2 * N + 1, 2, h)
	D2 = DM[:, :, 1]

	#	Reflect and add values appropriatey

	return (D2[-N-1:][:,-N-1:] + 
		hstack((zeros((N + 1, 1)), fliplr(D2[-N - 1:][:, 0:N]))))

def D1_x_e(N,h):
	"""
	Compute the (1/x)(d/dx) matrix for even functions.
	"""
	
	x = points(N, h)
	k = arange(1, N + 1)

	DM =  sincmats.deriv_mats(2 * N + 1, 1, h)
	D1 = DM[:, :, 0]

	
	#	This row needs to be specially constructed
	
	zero_row = hstack(([-pi**2 / (3 * h**2)], -4 * ((-1)**k) / ((h**2) * k**2)))
		
	#	Reflect and add values appropriatey
	
	return vstack(( zero_row, ( dot(diag(1.0 / x[-N:]), D1[-N:][:,-N-1:]) + 
		hstack(( zeros((N, 1)), fliplr(dot(diag(1.0 / x[-N:] ), D1[-N:][:, 0:N])) )) ) ))

def D1_eo(N,h):
	"""
	Compute the d1 matrix for mapping even functions to odd functions.
	"""
	
	DM = sincmats.deriv_mats(2 * N + 1, 1, h)
	D1 = DM[:, :, 0]

	#	Reflect and add values appropriatey

	return (D1[-N - 1:][:, -N - 1:] + 
		hstack((zeros((N + 1, 1)), fliplr(D1[-N - 1:][:, 0:N]))))
	
	
def IN1_oe(N,h):
	"""
	Compute the matrix of integration for mapping odd functions to even functions.
	"""
	
	IN1 = sincmats.quad_mat(2*N+1,h)

	#	Reflect and subtract values appropriatey

	return (IN1[-N - 1:][:,-N - 1:] - 
		hstack((zeros((N + 1, 1)), fliplr(IN1[-N - 1:][:, 0:N]))))
		
def sincinterp_e(x, u, xx):
	"""
	Interpolate at points xx for a function defined
	by the even vector u at the collocation points x.
	xx is assumed to lie inside [0,x(end)].
	"""
	
	h = x[1]-x[0]
	uu = zeros(xx.shape[0])
	
	for j in range(0, u.shape[0]):
		uu = uu + u[j]*sinc( (xx- x[j])/h)
		
	for j in range(1, u.shape[0]):
		uu = uu + u[j]*sinc( (xx + x[j])/h)
		
	return uu


