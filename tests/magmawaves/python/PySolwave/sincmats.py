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
Matices and collocation points for using the sinc spectral method.
Based on the MATLAB scripts DMSUITE by J.A.C. Weiderman.

Gideon Simpson
simpson@math.toronto.edu
"""


from numpy import linspace, imag, zeros, arange, pi, exp, concatenate
from scipy.linalg import toeplitz
from scipy.special import sici

def points(N, h):
	"""
	Compute the collocation points.
	"""

	return h*linspace(-(N - 1) / 2.0, (N - 1) / 2.0, N)

def deriv_mats(N,M,h):
	"""
	Compute the differentiation matrices.
	"""
	
	k = arange(1, N)
	t = pi*k
	sigma = zeros(N - 1)
	col = zeros(N)
	row = zeros(N)
	DM = zeros((N, N, M))
	
	for ell in range(1, M + 1):
		sigma = (-ell * sigma + imag(exp(1j * t) * (1j**ell))) / t
		col = (pi / h)**ell * concatenate( ([imag(1j**(ell + 1)) / (ell + 1)], sigma) )
		row = (-1)**ell*col
		row[0]=col[0]
		DM[:,:,ell-1] = toeplitz(col, row)
	
	return DM

def quad_mat(N, h):
	"""
	Compute the integration matrix.
	"""
	
	k = linspace(-(N - 1) / 2.0, (N - 1) / 2.0, N)	

	sigma = sici(pi * (k - k[0]))[0]
	col = h / 2.0 + (h / pi) * sigma
	row = h / 2.0 - (h / pi) * sigma
	
	return toeplitz(col, row)
