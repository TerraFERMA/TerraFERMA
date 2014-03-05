from dolfin import *
from solitarywave import *
from waveerrors import *


sw = SolitaryWave(5,3,0,2,400)

mesh = UnitSquare(64,64)
V = FunctionSpace(mesh,"CG",2)
f = ProjectSolitaryWave(V,[0.5,0.5],64.,sw)
plot(f)
interactive()

err = WaveError(f,sw,[0.5,0.5],64.)
print err.min_grad_z(0.)
