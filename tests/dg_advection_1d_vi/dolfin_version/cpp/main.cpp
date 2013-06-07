
#include <dolfin.h>

#include "advection.h"

using namespace dolfin;

class InitialCondition : public Expression
{
public:

  InitialCondition(boost::shared_ptr< MeshFunction<std::size_t> > cellids) : Expression(), cellids(cellids) {}

  void eval(Array<double>& values, const Array<double>& x,
            const ufc::cell& ufc_cell) const
  {
    const std::size_t cellid = (*cellids)[ufc_cell.index];
    if (cellid == 1)
      values[0] = 1.0;
    else
      values[0] = 0.0;
  }

private:

  boost::shared_ptr< MeshFunction<std::size_t> > cellids;

};


int main(int argc, char *argv[])
{
  Mesh mesh("interval.xml");
  boost::shared_ptr< MeshFunction<std::size_t> > cellids(new MeshFunction< std::size_t >(mesh, "interval_physical_region.xml"));
  boost::shared_ptr< MeshFunction<std::size_t> > facetids(new MeshFunction< std::size_t >(mesh, "interval_facet_region.xml"));

  // Create velocity FunctionSpace
  advection::FunctionSpace f_f(mesh);

  std::vector<double> vval(1, 4.0);
  Constant v(vval);
  Constant dt(0.01);
  Constant theta(0.5);
  Constant fin(1.0);

  Function f_i(f_f);
  Function f_n(f_f);
  InitialCondition f_0(cellids);
  f_i.interpolate(f_0);
  f_n.interpolate(f_0);

  advection::BilinearForm a(f_f, f_f);
  a.v = v; a.dt = dt; a.theta = theta;
  advection::LinearForm L(f_f);
  L.v = v; L.dt = dt; L.theta = theta; L.f_n = f_n; L.fin = fin;

  a.set_cell_domains(cellids);
  a.set_exterior_facet_domains(facetids);
  a.set_interior_facet_domains(facetids);
  L.set_cell_domains(cellids);
  L.set_exterior_facet_domains(facetids);
  L.set_interior_facet_domains(facetids);

  PETScMatrix A;
  PETScVector b;
  set_log_level(DEBUG);
  SystemAssembler assembler(a, L);
  assembler.assemble(A,b);
  assembler.reset_sparsity = false;

  double t = 0.0;
  double T = 1.0;

  File pvdfile("advection.pvd");
  pvdfile << std::make_pair<const Function*, double>(&f_i, t);

  while (t < T)
  {
    t += double(dt);
    solve(A, *f_i.vector(), b);
    *f_n.vector() = *f_i.vector();
    assembler.assemble(b);
    pvdfile << std::make_pair<const Function*, double>(&f_i, t);
  }

}
