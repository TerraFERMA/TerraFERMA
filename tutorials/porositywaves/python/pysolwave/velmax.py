# -*- coding: utf-8 -*-
"""
velmax.py:
    simple script to calculate maximum melt velocity at peak of a solitary wave
    for arbitrary c,m,n and d using Gideon's sinc collocation methods for 
    calculating solitary wave properties and derivatives.
    
    Here we calculate the maximum melt velocity at the peak of the wave of speed 
    c (amplitude A) as
    
    w = A^(n-1)*( 1 - dP/dz)
    
    given that in a solitary wave \phi^m P = - c d\phi/dz or
    
    w = A^(n-1)* ( 1 + c A^(-m) f'')
Created on Tue Feb  5 21:48:05 2013

@author: mspieg
"""

from .solitarywave import SolitaryWave
from .sinc_eo import D2_e
from numpy import *

import pylab as pl

def velmax(swave):
    h = diff(swave.r)[0]
    A = swave.f[0]
    n = swave.n
    c = swave.c
    m = swave.m
    d2fA = dot(D2_e(swave.N,h),swave.f-1.)
    w = A**(n-1)*( 1. + c*A**(-m)*d2fA[0])
    return A,w
    
def plotvelmax(c,w,n,m):
    pl.plot(c,w[0,:],'r')
    pl.hold(True)
    pl.plot(c,w[1,:],'b')
    pl.plot(c,w[2,:],'g')
    pl.plot(c,c,'k--')
    pl.legend(('d=1','d=2','d=3','c'),loc='best')
    pl.xlabel('wave speed c')
    pl.ylabel('v_max')
    pl.title("Maximum melt speed: n= {0:d}, m={1:d}".format(n,m))
    pl.grid()
    pl.hold(False)

def plotampmax(c,A,n,m):
    pl.plot(c,A[0,:],'r')
    pl.hold(True)
    pl.plot(c,A[1,:],'b')
    pl.plot(c,A[2,:],'g')
    pl.legend(('d=1','d=2','d=3'),loc='best')
    pl.xlabel('wave speed c')
    pl.ylabel('wave Amplitude')
    pl.title("Maximum Amplitude: n= {0:d}, m={1:d}".format(n,m))
    pl.grid()
    pl.hold(False)
    
# let's calculate w(c) for various waves
N = 300
m = 1
n = 2
nc = 20
c = linspace(1.1*n,10,nc)
Aar = zeros((3,nc))
war = zeros((3,nc))

for d in range(1,4):
    print('d = ',d)
    print('c     A   w')
    for i in range(len(c)):
        sw = SolitaryWave(c[i],n,m,d,N)
        A,w = velmax(sw)
        Aar[d-1,i] = A
        war[d-1,i] = w
        print(c[i],A,w)
        
pl.figure()
plotvelmax(c,war,n,m)
pl.figure()
plotampmax(c,Aar,n,m)
pl.show()

    
## quick test script
#d= 2
#m = 1
#n = 2
#N= 400
#
#c = 2.2
#sw = SolitaryWave(c,n,m,d,N)
#A,w = velmax(sw)
#
#h = diff(sw.r)[0]
#D2 = D2_e(N,h)
#d2f = dot(D2,sw.f-1)
#
#print 'c=',c,' n=',n,' m=',m,' d=',d,' A=',A,' w= ',w
#pl.plot(sw.r,sw.f,'r')
#pl.hold
#pl.plot(sw.r,d2f,'b')
#pl.legend(('f','D2f'),loc='best')
#pl.show()


