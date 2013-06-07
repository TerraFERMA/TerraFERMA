
from dolfin import *

class InitialCondition(Expression):
    def __init__(self, cellids):
        self.cellids = cellids
    def eval_cell(self, values, x, ufc_cell):
        cellid = self.cellids[ufc_cell.index]
        if cellid == 1:
          values[0] = 1.0
        else:
          values[0] = 0.0

mesh = Mesh("interval.xml")
cellids = MeshFunction("size_t", mesh, "interval_physical_region.xml")
facetids = MeshFunction("size_t", mesh, "interval_facet_region.xml")

f_f = FunctionSpace(mesh, "DG", 1)
v = Constant((1.0,))

dt = Constant(0.01)
dtr = dt('+')
fin = Constant(1.0)
theta = Constant(0.5)
thetar = theta('+')

f_i = Function(f_f)
f_n = Function(f_f)

f_0 = InitialCondition(cellids)
f_i.interpolate(f_0)
f_n.interpolate(f_0)

f_t = TestFunction(f_f)
f_a = TrialFunction(f_f)

n = f_f.cell().n
vin = (dot(v, n) + abs(dot(v, n)))/2.0
vout = (dot(v, n) - abs(dot(v, n)))/2.0
vn = dot(v,n)

r_m = f_t*f_a*dx(1) + f_t*f_a*dx(2) - f_t*f_n*dx(1) - f_t*f_n*dx(2)

r_a = - dt*theta*dot(grad(f_t), v*f_a)*dx(1) - dt*theta*dot(grad(f_t), v*f_a)*dx(2) \
      - dt*(1.0-theta)*dot(grad(f_t), v*f_n)*dx(1) - dt*(1.0-theta)*dot(grad(f_t), v*f_n)*dx(2)

#r_m = f_t*f_a*dx - f_t*f_n*dx
#
#r_a = - dt*theta*dot(grad(f_t), v*f_a)*dx - dt*(1.0-theta)*dot(grad(f_t), v*f_n)*dx

r_fS =  dtr*thetar*(dot(vin('+')*f_a('+') - vin('-')*f_a('-'), jump(f_t))*dS) + \
        dtr*(1.0-thetar)*(dot(vin('+')*f_n('+') - vin('-')*f_n('-'), jump(f_t))*dS)

r_fs1 =  dt*dot(vn*f_t, fin)*ds(1)
r_fs2 =  dt*theta*dot(vn*f_t, f_a)*ds(2) + dt*(1.0-theta)*dot(vn*f_t, f_n)*ds(2)
r_fs = r_fs1 + r_fs2

#r_fs =  dt*dot(vout*f_t, fin)*ds + dt*theta*dot(vin*f_t, f_a)*ds + dt*(1.0-theta)*dot(vin*f_t, f_n)*ds

r = r_m + r_a + r_fS + r_fs

a = Form(lhs(r)) # seem to need to use this copy constructor to get the option of
L = Form(rhs(r)) # adding the domain info
a.set_cell_domains(cellids)
L.set_cell_domains(cellids)
a.set_interior_facet_domains(facetids)
L.set_interior_facet_domains(facetids)
a.set_exterior_facet_domains(facetids)
L.set_exterior_facet_domains(facetids)
# the above (attaching of domain ids) doesn't seem to work

set_log_level(DEBUG)
assembler = SystemAssembler(a, L)

A = PETScMatrix()
b = PETScVector()

assembler.assemble(A, b)
assembler.reset_sparsity = False

pvdfile = File("advection.pvd")

t = 0.0
T = 1.0
pvdfile << (f_i, t)

while (t < T):
  t += float(dt)
  solve(A, f_i.vector(), b)
  f_n.vector()[:] = f_i.vector()
  assembler.assemble(b)
  pvdfile << (f_i, t)



