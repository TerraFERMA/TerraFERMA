"""python functions for calculating solitary waves and errors in dolfin functions
"""

__author__ = "Marc Spiegelman (mspieg@ldeo.columbia.edu)"
__date__ = "15 Jun 2010 11:39:55"
__copyright__ = "Copyright (C) 2010 Marc Spiegelman"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
from numpy import *
from scipy.optimize import fmin_bfgs, fsolve, brent, brentq
import sys
from ttictoc import tic,toc

from .magmasinc import solwave_mck, solwave_con, solwave_m1, solwave_gen
from .solitarywave import *

def radial_expression(wd,d):
    """ returns Expression for radial function R= || x - x_0|| for wd dimensional solitary wave in d dimensional space
    """

    # set up "radial function" in wd dimensions 
    if wd == 1: 
        str = "sqrt(pow((x[%d]-x%d),2))" % (d-1,d-1)
    elif wd == 2:
        str = "sqrt(pow((x[0]-x0),2) + pow((x[%d]-x%d),2))" % (d-1,d-1)
    elif wd == 3:
        str = "sqrt(pow((x[0]-x0),2) + pow((x[1]-x1),2) + pow((x[2]-x2),2))"


    R=Expression(str,x0=0.,x1=0.,x2=0.,degree=1)
    return R


def ProjectSolitaryWave(V,x0,h_on_delta,solitarywave):
    """ Calculate the projection of a Solitary Wave onto a function space V

    x0: origin of solitary wave on domain
    h_on_delta: scale length in compaction lengths
    solitarywave:  SolitaryWave profile object (from solitarywave.py)
    """


    # initialize the radial Expression function
    dim = V.mesh().geometry().dim()
    R = radial_expression(solitarywave.d,dim)
    R.x0 = x0[0]
    if dim >= 2:
        R.x1 = x0[1]
    if dim == 3:
        R.x2 = x0[2]

    # project r onto Function f
    f = project(R,V)

    # scale R by h_on_delta
    fx = f.vector()
    fx *= h_on_delta
    
    #get numpy array for radial function and interpolate using sinc
    rmesh = array([f for f in fx]) # get numpy array
    phi = solitarywave.eval(rmesh)
    
    # reload phi onto function
    for i in range(fx.size()):
        fx[i] = phi[i]

    return f
            
    

class WaveError:
    """Class for calculating L2 errors for solitary wave problems"""

    
    def __init__(self, f,solitarywave,x0,h_on_delta):
        """Initialize Error Object

        f: the Discrete Finite Element solution for porosity
        solitarywave: analtyic solitary wave object (from class solitarywave)
        x0: origin of solitary wave on domain
        h_on_delta: scale length in compaction lengths
        
        """

        # the Discrete Finite Element solution for porosity
        self.f = f
        # the analytic solitarywave profile object
        self.solitarywave = solitarywave

        # Calculate and initialize radial function Expression, on P2 CG elements
        dim = f.function_space().mesh().geometry().dim()
        self.dim = dim
        self.x0 = x0
        self.R = radial_expression(solitarywave.d,dim)
        self.setRx(self.x0)

        # set scaling
        self.h_on_delta = h_on_delta

        #Initialize Quadrature Function for fc
        #mesh = f.function_space().mesh()
        #QE = FunctionSpace(mesh,"Quadrature",4)
        #QE = FunctionSpace(mesh,"CG",2) # try P2 for the moment
        self.fc = Function(f.function_space())
        
        #Initialize Form for Error Functional
        err = (self.f-self.fc)
        self.M = err*err*dx
        M = self.f*self.f*dx

        #calculate norm of function for relative errors
        #L2_f = assemble(M, mesh = f.function_space().mesh())
        L2_f = assemble(M)
        self.L2_f = sqrt(L2_f)

        #calculate inner product 2<f-fc,Grad(fc)>
        self.M0 = 2*err*self.fc.dx(0)*dx
        self.M1 = 2*err*self.fc.dx(1)*dx
        if dim == 3:
            self.M2 = 2*err*self.fc.dx(2)*dx

    def setRx(self,x):
        """ set the origin of the radial function to x
        """
        self.R.x0 = x[0]
        if self.dim >= 2:
            self.R.x1 = x[1]
        if self.dim == 3:
            self.R.x2 = x[2]

    def set_fc(self,delta):
        """ interpolate the exact solitary wave with offset x0 + delta onto a function in Function Space V
        """
        # Interpolate R onto the function
        self.setRx(self.x0+delta)
        self.fc.interpolate(self.R)

        # scale R by h_on_delta
        fx = self.fc.vector()
        fx *= self.h_on_delta

        #get numpy array for radial function and interpolate using sinc
        rmesh = array([f for f in fx]) # get numpy array
        phi = self.solitarywave.eval(rmesh)

        # reload phi onto function
        for i in range(fx.size()):
            fx[i] = phi[i]
      
  
    def L2(self,delta):
        """Calculate the L2 norm of the error ||f - fc(delta)||
        """

        # set the solitary wave fc with offset delta
        self.set_fc(delta)
        L2_error = assemble(self.M)

        return sqrt(L2_error)

    def gradL2(self,delta):
        """returns the gradient of the error functional with respect to delta evaluated at delta
        """
        # set the solitary wave fc with offset delta
        self.set_fc(delta)
        mesh = self.f.function_space().mesh()
        L2_0 = assemble(self.M0)
        L2_1 = assemble(self.M1)
        if self.dim < 3:
            return array([L2_0, L2_1])
        elif self.dim == 3:
            L2_2 = assemble(self.M2)
            return array([L2_0, L2_1, L2_2])

    def gradL2_z(self,delta_z):
        """just returns the z component of the gradient of the error functional with respect to delta evaluated at delta
        assumes delta_x, delta_y =0
        """
        # set the solitary wave fc with offset delta_z
        delta_1d = zeros(self.dim)
        delta_1d[-1] = delta_z;
        
        self.set_fc(delta_1d)
        mesh = self.f.function_space().mesh()
        if self.dim ==2:
            L2_z = assemble(self.M1)
        elif self.dim == 3:
            L2_z = assemble(self.M2)
            
        return L2_z
        

    def min(self,delta_0): 
        """uses a bfgs optimizer to minimize L2(delta)
        """
        tic()
        self.fopt = fmin_bfgs(self.L2,delta_0,fprime=self.gradL2,disp=1,retall=1,full_output=1)
        elapsed_time = toc()

        delta_min = self.fopt[0]
        L2_min = self.fopt[1]

        return [delta_min, L2_min, elapsed_time]

    def min_grad(self,delta_0):
        """uses fsolve (without jacobians) to find the roots of Grad E =0, then finds minimum
        """

        tic()
        self.fopt = fsolve(self.gradL2,delta_0,full_output = 1)
        elapsed_time = toc()

        delta_min =self.fopt[0]
        L2_min = self.L2(delta_min)

        return [delta_min, L2_min, elapsed_time]

    def min_grad_z(self,delta_0):
        """uses fbrentq to just find the root for dE/d_delta_z = 0 (assumes delta_x = 0 )
        """

        tic()
        self.fopt = brentq(self.gradL2_z,delta_0-.2,delta_0+.2,full_output = 1)
        elapsed_time = toc()

        delta_z_min =self.fopt[0]
        delta_min = zeros(self.dim)
        delta_min[-1] = delta_z_min
        L2_min = self.L2(delta_min)

        return [delta_min, L2_min, elapsed_time]

