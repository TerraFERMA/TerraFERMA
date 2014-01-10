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

"""python class  for calculating solitary waves profiles
"""

__author__ = "Marc Spiegelman (mspieg@ldeo.columbia.edu)"
__date__ = "15 Oct 2011 22:01:49"
__copyright__ = "Copyright (C) 2010 Marc Spiegelman"
__license__  = "GNU LGPL Version 2.1"

from numpy import *
from scipy.interpolate import interp1d
import sys

from magmasinc import solwave_mck, solwave_con, solwave_m1, solwave_gen
from sinc_eo import sincinterp_e

class SolitaryWave:
    """ class for calculating and evaluating solitary wave profiles from Simpson
    magmasinc library
    """

    def __init__(self, c,n,m,d,N):
        """initialize radial solitary wave profile fc(r) using PySolwave Routines

        c: wavespeed
        n: permeability exponent
        m: bulk viscosity exponent
        d: wave dimension
        N: number of collocation points
        """

        # calculate profiles checking for special cases
        if m == 1:
            if n == 2:
                r, f = solwave_con(c, d, N)
            else:
                r, f = solwave_m1(c,n, d, N)
        else:
            if n == 3:
                r, f = solwave_mck(c, d, N)
            else:
                r, f = solwave_gen(c, n, m, d, N)
                
        self.r = r
        self.f = f
        self.finterp = interp1d(r,f, kind='cubic', fill_value=1.0)

    def eval(self,r):
        """ use sincinterp_e to return porosity at radius r"""

        f = sincinterp_e(self.r, self.f - 1.0, r) + 1.0
        return f


    def interp(self,r):
        """ use scipy interp1d cubic interpolation to return values"""
        return self.finterp(r)

    

        
        

                
