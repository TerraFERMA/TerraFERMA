
#ifndef __DOLFIN_BOOST_TYPES_H
#define __DOLFIN_BOOST_TYPES_H

#include <dolfin.h>
#include "petscsnes.h"

namespace buckettools {

  // predeclarations
  class SystemBucket;
  typedef boost::shared_ptr< SystemBucket > SystemBucket_ptr;
  class FunctionBucket;
  typedef boost::shared_ptr< FunctionBucket > FunctionBucket_ptr;
  class SolverBucket;
  typedef boost::shared_ptr< SolverBucket > SolverBucket_ptr;

  // Define boost shared_ptr types for lots of useful things from DOLFIN
  typedef boost::shared_ptr< dolfin::GenericFunction >              GenericFunction_ptr;
  typedef boost::shared_ptr< dolfin::Expression >                   Expression_ptr;
  typedef boost::shared_ptr< dolfin::Mesh >                         Mesh_ptr;
  typedef boost::shared_ptr< dolfin::MeshFunction< dolfin::uint > > MeshFunction_uint_ptr;
  typedef boost::shared_ptr< dolfin::FunctionSpace >                FunctionSpace_ptr;
  typedef boost::shared_ptr< dolfin::Function >                     Function_ptr;
  typedef boost::shared_ptr< dolfin::BoundaryCondition >            BoundaryCondition_ptr;
  typedef boost::shared_ptr< dolfin::DirichletBC >                  DirichletBC_ptr;
  typedef boost::shared_ptr< dolfin::Form >                         Form_ptr;
  typedef boost::shared_ptr< dolfin::PETScMatrix >                  PETScMatrix_ptr;
  typedef boost::shared_ptr< dolfin::PETScVector >                  PETScVector_ptr;

  typedef boost::shared_ptr< KSP > KSP_ptr;
  typedef boost::shared_ptr< PC >  PC_ptr;
  typedef boost::shared_ptr< IS >  IS_ptr;
  typedef boost::shared_ptr< Mat > Mat_ptr;
  typedef boost::shared_ptr< Vec > Vec_ptr;

  // Define iterator types for things accessed in the bucket maps (defined below)
  typedef std::map< std::string, SystemBucket_ptr >::iterator           SystemBucket_it;
  typedef std::map< std::string, SystemBucket_ptr >::const_iterator     SystemBucket_const_it;
  typedef std::map< int, SystemBucket_ptr >::iterator                   int_SystemBucket_it;
  typedef std::map< int, SystemBucket_ptr >::const_iterator             int_SystemBucket_const_it;
  typedef std::map< std::string, Mesh_ptr >::iterator                   Mesh_it;
  typedef std::map< std::string, Mesh_ptr >::const_iterator             Mesh_const_it;
  typedef std::map< std::string, std::string >::iterator                string_it;
  typedef std::map< std::string, std::string >::const_iterator          string_const_it;
//  typedef std::map< std::string, GenericDetectors_ptr >::iterator       GenericDetectors_it;
//  typedef std::map< std::string, GenericDetectors_ptr >::const_iterator GenericDetectors_const_it;
  
  // Define iterator types for things accessed in the system maps (defined below)
  typedef std::map< std::string, FunctionSpace_ptr >::iterator           FunctionSpace_it;
  typedef std::map< std::string, FunctionSpace_ptr >::const_iterator     FunctionSpace_const_it;
  typedef std::map< std::string, FunctionBucket_ptr >::iterator          FunctionBucket_it;
  typedef std::map< std::string, FunctionBucket_ptr >::const_iterator    FunctionBucket_const_it;
  typedef std::map< std::string, SolverBucket_ptr >::iterator            SolverBucket_it;
  typedef std::map< std::string, SolverBucket_ptr >::const_iterator      SolverBucket_const_it;
  typedef std::map< std::string, GenericFunction_ptr >::iterator         GenericFunction_it;
  typedef std::map< std::string, GenericFunction_ptr >::const_iterator   GenericFunction_const_it;
  typedef std::map< std::string, Function_ptr >::iterator                Function_it;
  typedef std::map< std::string, Function_ptr >::const_iterator          Function_const_it;
  typedef std::map< std::string, Expression_ptr >::iterator              Expression_it;
  typedef std::map< std::string, Expression_ptr >::const_iterator        Expression_const_it;
  typedef std::map< std::string, DirichletBC_ptr >::iterator             DirichletBC_it;
  typedef std::map< std::string, DirichletBC_ptr >::const_iterator       DirichletBC_const_it;
  typedef std::map< std::string, BoundaryCondition_ptr >::iterator       BoundaryCondition_it;
  typedef std::map< std::string, BoundaryCondition_ptr >::const_iterator BoundaryCondition_const_it;
  typedef std::map< uint, Expression_ptr >::iterator                     uint_Expression_it;
  typedef std::map< uint, Expression_ptr >::const_iterator               uint_Expression_const_it;
  typedef std::map< int, FunctionSpace_ptr >::iterator                   int_FunctionSpace_it;
  typedef std::map< int, FunctionSpace_ptr >::const_iterator             int_FunctionSpace_const_it;

  // Define iterator types for things accessed in the function maps (defined below)
  typedef std::map< std::string, Form_ptr >::iterator       Form_it;
  typedef std::map< std::string, Form_ptr >::const_iterator Form_const_it;
  
}

#endif
