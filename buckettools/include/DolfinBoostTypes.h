
#ifndef __DOLFIN_BOOST_TYPES_H
#define __DOLFIN_BOOST_TYPES_H

#include <dolfin.h>

namespace buckettools {

  // Define boost shared_ptr types for lots of useful things
  typedef boost::shared_ptr< dolfin::GenericFunction >              GenericFunction_ptr;
  typedef boost::shared_ptr< dolfin::Expression >                   Expression_ptr;
  typedef boost::shared_ptr< dolfin::Mesh >                         Mesh_ptr;
  typedef boost::shared_ptr< dolfin::MeshFunction< dolfin::uint > > MeshFunction_uint_ptr;
  typedef boost::shared_ptr< dolfin::FunctionSpace >                FunctionSpace_ptr;
  typedef boost::shared_ptr< dolfin::Function >                     Function_ptr;
  typedef boost::shared_ptr< dolfin::DirichletBC >                  DirichletBC_ptr;
  typedef boost::shared_ptr< dolfin::Form >                         Form_ptr;

}

#endif
