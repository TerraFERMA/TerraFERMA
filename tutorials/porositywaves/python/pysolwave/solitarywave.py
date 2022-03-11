"""python class  for calculating solitary waves profiles
"""

__author__ = "Marc Spiegelman (mspieg@ldeo.columbia.edu)"
__date__ = "15 Oct 2011 22:01:49"
__copyright__ = "Copyright (C) 2010 Marc Spiegelman"
__license__  = "GNU LGPL Version 2.1"

from numpy import *
from scipy.interpolate import interp1d
import sys

from .magmasinc import solwave_mck, solwave_con, solwave_m1, solwave_gen
from .sinc_eo import sincinterp_e

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
        # initialize wave parameters
        self.c = c
        self.n = n
        self.m = m
        self.d = d
        self.N = N
        
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

    

        
        

                
