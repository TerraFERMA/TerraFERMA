from pysolwave.solitarywave import *
import numpy as np
import matplotlib.pyplot as pl

c = 5 # wave speed
n = 3 # permeability exponent K=\phi^n
m = 1 # bulk porosity exponent \zeta = \phi^{-m}
N = 150 # number of collocation points for sinc-method
dim = [1,2,3] # list of dimensions
swa = [] #  list of solitary wave objects

# calculate and store solitary wave objects
for d in dim:
    sw = SolitaryWave(c,n,m,d,N)
    swa.append(sw)
    
# and plot them out
pl.figure()
for i in range(len(dim)):
    pl.plot(swa[i].r,swa[i].f,'+-',label='d={0}'.format(dim[i]))
    pl.hold(True)

pl.xlabel('Radial distance r')
pl.ylabel('porosity f')
pl.ylim((0.5, 4.))
pl.legend(loc='best')
pl.title('Solitary Wave profiles c={0}, n={1}, m={2}'.format(c,n,m))
pl.grid()
pl.savefig('SolitaryWavesProfiles.pdf')
pl.show()
