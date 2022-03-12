#!/usr/bin/env python
# script to extract phase and shape errors from 2-D magmawave benchmarks
# Marc Spiegelman

__author__ = "Marc Spiegelman (mspieg@ldeo.columbia.edu)"
__date__ = " 14 Mar 2012 15:51:51"
__copyright__ = "Copyright (C) 2010 Marc Spiegelman"
__license__  = "GNU LGPL Version 3.1"

from dolfin import *
import os, sys, subprocess
from glob import glob
from pickle import dump

# note:  Current PySolwave must be in python path

from .solitarywave import *
from .waveerrors import *

# and need libspud
import libspud

def tfmlgeterrors(name):
    # get arguments from .tfml file
    tfmlfile = name+'.tfml'
    libspud.load_options(tfmlfile)
    
    
    # get the mesh parameters and reconstruct the mesh
    meshtype = libspud.get_option("/geometry/mesh/source[0]/name")
    if meshtype == 'UnitSquare':
        number_cells = libspud.get_option("/geometry/mesh[0]/source[0]/number_cells")
        diagonal = libspud.get_option("/geometry/mesh[0]/source[0]/diagonal")
        mesh = UnitSquare(number_cells[0],number_cells[1],diagonal)
    elif meshtype == 'UnitCube':
        number_cells = libspud.get_option("/geometry/mesh[0]/source[0]/number_cells")
        mesh = UnitCube(number_cells[0],number_cells[1],number_cells[2])
    else:
       error("Error: unknown mesh type "+meshtype)

    # get the mesh dimension 
    dim = mesh.topology().dim()
    
    # get and set  the mixed function space
    p_family = libspud.get_option("/system[0]/field[0]/type/rank/element/family")
    p_degree = libspud.get_option("/system[0]/field[0]/type/rank/element/degree")
    f_family = libspud.get_option("/system[0]/field[1]/type/rank/element/family")
    f_degree = libspud.get_option("/system[0]/field[1]/type/rank/element/degree")
    Vp = FunctionSpace(mesh,p_family,p_degree)
    Vf = FunctionSpace(mesh,f_family,f_degree)
    ME = Vp*Vf
    
    # get the solitarywave profile
    solwave_code = libspud.get_option("/system::magma/field::Porosity/type::Function/rank::Scalar/initial_condition/python")
    exec(solwave_code)
    # print wave parameters
    print('Solitary Wave parameters: c =',swave.c,', n = ',swave.n,', m = ',swave.m,', d = ',swave.d)
    print('h_on_delta =', h_on_delta)
   
    # calculate dh_nodes_per_compaction length
    dh = mesh.hmin() # cell diameter (circum-circle)
    dh_node = dh/sqrt(dim)/2.
    dh_node_delta = dh_node*h_on_delta

    # cd to build directory to loop over checkpoint files
    # os.chdir('build')

    # get and sort xml files by time
    tfml_files = glob(name+'*checkpoint*.tfml')
    xml_files = glob(name+'*.xml')
    dtype = [ ('time', float), ('name','S28')]
    values = []
    for i in range(len(tfml_files)):
        libspud.load_options(tfml_files[i])
        time = libspud.get_option("/timestepping/current_time")
        dt = libspud.get_option("/timestepping/timestep/coefficient/type/rank/value/constant")
        values = values + [(time, xml_files[i])]

    a = array(values,dtype)
    t_file=sort(a,order='time')
         
    
    set_log_active(False) # turn off dolfin messages

    # initialize output arrays
    t = []
    delta = []
    err_min = []
    L2_f = []
    rel_error = []
    c_rel_error = []
    elapsed_time = []
    
    # loop over files and check errors for each one

    delta_0 = zeros(dim)
    for i in range(len(t_file)):
        
        # get time and filename
        t.append(t_file[i][0])
        name = t_file[i][1]
        
        # load function and extract porosity
        u = Function(ME,name)
        p, f = u.split(deepcopy = True)
        
        
        #initialize Error Object
        err = WaveError(f,swave,x0,h_on_delta)
        
        #minimize error using fsolve
        #out = err.min_grad(delta_0)
        out = err.min_grad_z(delta_0[-1])

        # collect output variables
        delta_i = out[0]
        delta.append(delta_i)
        err_min.append(out[1])
        elapsed_time = out[2]
        #print "delta= ", delta_i

        L2_f.append(err.L2_f)
        rel_error.append(err_min[-1]/L2_f[-1])
        c_rel_error.append(sign(delta_i[-1])*sqrt(dot(delta_i,delta_i))/(c*t[-1]))

        #create string and dump to screen and file
        if dim == 2:
            str = "%g\t [ %g %g ] \t %g\t %g\t %g\t %g\t%g" % (t[-1],delta_i[0],delta_i[1],err_min[-1],L2_f[-1],rel_error[-1],c_rel_error[-1],elapsed_time)
        elif dim == 3:
            str = "%g\t [ %g %g %g ] \t %g\t %g\t %g\t %g\t%g" % (t[-1],delta_i[0],delta_i[1],delta_i[2],err_min[-1],L2_f[-1],rel_error[-1],c_rel_error[-1],elapsed_time)

                
        disp(str)
        delta_0 = delta_i

    # set up dictionary and pickle
    list = [['c', swave.c],['d', swave.d],['n',swave.n],['m',swave.m],
            ['h_on_delta',h_on_delta],['dh_on_delta', dh_node_delta],['dt',dt],
            [ 'time', t], ['delta', delta], ['min_err', err_min], ['L2_f', L2_f], ['rel_error', rel_error], ['c_rel_error',c_rel_error]]
    d = dict(list)
    f = file('ErrorFile.pkl','w+')
    dump(d,f)

