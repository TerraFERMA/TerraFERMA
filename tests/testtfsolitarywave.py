# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 22:00:36 2015

@author: mspieg
"""

import numpy as np
from glob import glob
from pysolwave.tfsolitarywave import *
import pylab as pl

# load the tfml fie
tfname = glob("*.tfml")[0]
tf = TFSolitaryWave(tfname)


x0 = tf.x0
print "wave origin=",x0

r = tf.getr(x0)
f = tf.eval(x0)
print "x=",x0," r=",r, "f=", f

l = np.linspace(1,2,3)
x = np.outer(l,x0)
r = tf.getr(x)
f = tf.eval(x)

print "x=",x," r=",r, "f=", f

x[:,0] = x0[0]
r = tf.getr(x)
f = tf.eval(x)

print "x=",x," r=",r, "f=", f
pl.figure()
pl.plot(r,f,'bo',tf.swave.r,tf.swave.f,'r--')

checkpoint=glob("*.tfml")[1]
print "using checkpoint file ", checkpoint
errors = tf.geterrors(checkpoint)
print errors

