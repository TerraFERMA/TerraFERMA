from dolfin import *
from pysolwave.solitarywave import *
from pysolwave.waveerrors import *


sw = SolitaryWave(5,3,0,2,400)

mesh = UnitSquareMesh(32,32)
V = FunctionSpace(mesh,"CG",2)
f = ProjectSolitaryWave(V,[0.5,0.6],64.,sw)
plot(f)
interactive()

err = WaveError(f,sw,[0.5,0.5],64.)
print err.min_grad_z(0.)
