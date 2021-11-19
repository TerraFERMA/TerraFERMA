"""python class  for calculating solitary waves given input from TerraFERMA .tfml files
   uses libspud to read the tfml files and then populate
"""

__author__ = "Marc Spiegelman (mspieg@ldeo.columbia.edu)"
__date__ = "8 February 2015"
__copyright__ = "Copyright (C) 2010 Marc Spiegelman"
__license__  = "GNU LGPL Version 2.1"

import numpy as np
from .solitarywave import SolitaryWave
from .waveerrors import WaveError
import libspud
import dolfin as df

class TFSolitaryWave:
    """ class for calculating and evaluating solitary wave profiles with data from TerraFERMA .tfml input files
    """

    def __init__(self, tfml_file,system_name='magma',p_name='Pressure',f_name='Porosity',c_name='c',n_name='n',m_name='m',d_name='d',N_name='N',h_squared_name='h_squared',x0_name='x0'):
        """read the tfml_file and use libspud to populate the internal parameters

        c: wavespeed
        n: permeability exponent
        m: bulk viscosity exponent
        d: wave dimension
        N: number of collocation points
        x0: coordinate wave peak
        h_squared:  the size of the system in compaction lengths
                    squared (h/delta)**2
        """
        # initialize libspud and extract parameters
        libspud.clear_options()
        libspud.load_options(tfml_file)
        # get model dimension
        self.dim = libspud.get_option("/geometry/dimension")
        self.system_name = system_name
        # get solitary wave parameters
        path="/system::"+system_name+"/coefficient::"
        scalar_value="/type::Constant/rank::Scalar/value::WholeMesh/constant"
        vector_value="/type::Constant/rank::Vector/value::WholeMesh/constant::dim"
        c = libspud.get_option(path+c_name+scalar_value)
        n = int(libspud.get_option(path+n_name+scalar_value))
        m = int(libspud.get_option(path+m_name+scalar_value))
        d = float(libspud.get_option(path+d_name+scalar_value))
        N = int(libspud.get_option(path+N_name+scalar_value))
        self.h = np.sqrt(libspud.get_option(path+h_squared_name+scalar_value))
        self.x0 = np.array(libspud.get_option(path+x0_name+vector_value))
        self.swave = SolitaryWave(c,n,m,d,N)
        self.rmax = self.swave.r[-1]
        self.tfml_file = tfml_file
  
        # check that d <= dim
        assert (d <= self.dim)   
        
        # sort out appropriate index for calculating distance r
        if d == 1:
            self.index = [self.dim - 1]
        else: 
            self.index = list(range(0,int(d)))  
            
        # check that the origin point is the correct dimension
        assert (len(self.x0) == self.dim)
        
        
        #read in information for constructing Function space and dolfin objects
        # get the mesh parameters and reconstruct the mesh
        meshtype = libspud.get_option("/geometry/mesh/source[0]/name")
        if meshtype == 'UnitSquare':
            number_cells = libspud.get_option("/geometry/mesh[0]/source[0]/number_cells")
            diagonal = libspud.get_option("/geometry/mesh[0]/source[0]/diagonal")
            mesh = df.UnitSquareMesh(number_cells[0],number_cells[1],diagonal)
        elif meshtype == 'Rectangle':
            x0 = libspud.get_option("/geometry/mesh::Mesh/source::Rectangle/lower_left")
            x1 = libspud.get_option("/geometry/mesh::Mesh/source::Rectangle/upper_right")
            number_cells = libspud.get_option("/geometry/mesh::Mesh/source::Rectangle/number_cells")
            diagonal = libspud.get_option("/geometry/mesh[0]/source[0]/diagonal")
            mesh = df.RectangleMesh(x0[0],x0[1],x1[0],x1[1],number_cells[0],number_cells[1],diagonal)
        elif meshtype == 'UnitCube':
            number_cells = libspud.get_option("/geometry/mesh[0]/source[0]/number_cells")
            mesh = df.UnitCubeMesh(number_cells[0],number_cells[1],number_cells[2])
        elif meshtype == 'Box':
            x0 = libspud.get_option("/geometry/mesh::Mesh/source::Box/lower_back_left")
            x1 = libspud.get_option("/geometry/mesh::Mesh/source::Box/upper_front_right")
            number_cells = libspud.get_option("/geometry/mesh::Mesh/source::Box/number_cells")
            mesh = df.BoxMesh(x0[0],x0[1],x0[2],x1[0],x1[1],x1[2],number_cells[0],number_cells[1],number_cells[2])
        elif meshtype == 'UnitInterval':
            number_cells = libspud.get_option("/geometry/mesh::Mesh/source::UnitInterval/number_cells")
            mesh = df.UnitIntervalMesh(number_cells)
        elif meshtype == 'Interval':
            number_cells = libspud.get_option("/geometry/mesh::Mesh/source::Interval/number_cells")
            left = libspud.get_option("/geometry/mesh::Mesh/source::Interval/left")
            right = libspud.get_option("/geometry/mesh::Mesh/source::Interval/right")
            mesh = df.IntervalMesh(number_cells,left,right)
        elif meshtype == 'File':
            mesh_filename = libspud.get_option("/geometry/mesh::Mesh/source::File/file")
            print("tfml_file  = ",self.tfml_file, "mesh_filename=",mesh_filename)
            mesh = df.Mesh(mesh_filename)
        else:
            df.error("Error: unknown mesh type "+meshtype)
           
        #set the functionspace for n-d solitary waves
        path="/system::"+system_name+"/field::"
        p_family = libspud.get_option(path+p_name+"/type/rank/element/family")
        p_degree = libspud.get_option(path+p_name+"/type/rank/element/degree")
        f_family = libspud.get_option(path+f_name+"/type/rank/element/family")
        f_degree = libspud.get_option(path+f_name+"/type/rank/element/degree")        
        pe = df.FiniteElement(p_family, mesh.ufl_cell(), p_degree)
        ve = df.FiniteElement(f_family, mesh.ufl_cell(), f_degree)
        e = pe*ve
        self.functionspace = df.FunctionSpace(mesh, e)

        #work out the order of the fields
        for i in range(libspud.option_count("/system::"+system_name+"/field")):
          name = libspud.get_option("/system::"+system_name+"/field["+repr(i)+"]/name")
          if name == f_name:
            self.f_index = i
          if name == p_name:
            self.p_index = i
                
    def getr(self,x):
        """ return radial position with respect to wave origin x0
        for solitarywave of dimension d in overall dimension dim
        assuming x is a numpy array of points of dimension d"""
        
        # check that x is the right shaped numpy array
        
        if len(x.shape) == 1:             # 1-D array
            if self.dim == 1:
                npnts = x.shape[0]
            else:
                assert(x.shape == self.x0.shape)
                npnts = 1
        else:  #n-D array
            assert(x.shape[-1] == self.x0.shape[0])
            npnts = x.shape[0]
                    
        if npnts == 1:
            dx = x[self.index] - self.x0[self.index]
            r = self.h*np.sqrt(np.sum(dx*dx))
        else:
            dx = x[:,self.index] - self.x0[self.index]
            r = self.h*np.sqrt(np.sum(dx*dx,1))
        
        return r

    def eval(self,x):
        """ calculate position r given numpy array of x-coordinates
            and return appropriate solitary wave amplitude
        """
        

        # calculate radial coordinate r of x
        r = self.getr(x)

        # check if r is a scalar
        if np.isscalar(r):
            if r > self.rmax:
                f = 1.
            else:
                f = self.swave.interp(np.array([r]))
        else:
            f = np.ones(r.shape)
            omega = r <= self.rmax # r within collocation domain
            f[omega] = self.swave.finterp(r[omega])
        
        return f
        
    def geterrors(self,checkpoint_file):
        """
        reads checkpoint file and accompanying .xml solution file and
        extracts the time, and calculates shape and phase errors for the solitary waves
        """
        # this seems to be a dangerous thing but I'm going to clear
        # the option then load the new ones from the checkpoint file
        libspud.clear_options()
        libspud.load_options(checkpoint_file)
        time = libspud.get_option("/timestepping/current_time")
        dt = libspud.get_option("/timestepping/timestep/coefficient/type/rank/value/constant")
        print("t=",time," dt=", dt)
        
        #load xml file
        xml_file = checkpoint_file.replace("checkpoint",self.system_name)
        xml_file = xml_file.replace("tfml","xml")
        
        
        # load function and extract porosity
        u = df.Function(self.functionspace,xml_file)
        fields = u.split(deepcopy = True)
        f = fields[self.f_index]
        
        
        #initialize Error Object
        err = WaveError(f,self.swave,self.x0,self.h)
        
        #minimize error using fsolve
        #out = err.min_grad(delta_0)
        delta_0=np.zeros(self.dim)
        out = err.min_grad_z(delta_0[-1])
        
        delta = out[0]
        err_min = out[1]
        elapsed_time = out[2]
        
        L2_f= err.L2_f
        rel_error = err_min/L2_f
        c_rel_error = np.sign(delta)*np.sqrt(np.dot(delta,delta))/(self.swave.c*time)

        #print "delta=",delta, " err =",err_min, " elapsed_time=",elapsed_time
        #print "L2_f=",L2_f, "rel_err =",rel_error, " c_rel_error=",c_rel_error
        #print "L2_error=",err_min," delta=",delta," rel_error=",rel_error," c_rel_error=",c_rel_error
        return [time,dt,err_min,delta,rel_error,c_rel_error,elapsed_time]
        
        
