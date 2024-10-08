// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __DOLFIN_BOOST_TYPES_H
#define __DOLFIN_BOOST_TYPES_H

#include <dolfin.h>
#include "petscsnes.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

namespace buckettools {

  //*****************************************************************|************************************************************//
  // A collection of typedefs for std pointers to dolfin, petsc and bucket structures
  //*****************************************************************|************************************************************//

  //*****************************************************************|************************************************************//
  // std shared pointers to basic objects
  //*****************************************************************|************************************************************//

  typedef std::shared_ptr< int >    int_ptr;
  typedef std::shared_ptr< double > double_ptr;
  typedef std::shared_ptr< bool >   bool_ptr;

  //*****************************************************************|************************************************************//
  // std shared pointers to bucket objects
  //*****************************************************************|************************************************************//

  class Bucket;                                                      // predeclaration
  typedef std::shared_ptr< Bucket >           Bucket_ptr;
  class SystemBucket;                                                // predeclaration
  typedef std::shared_ptr< SystemBucket >     SystemBucket_ptr;
  class FunctionBucket;                                              // predeclaration
  typedef std::shared_ptr< FunctionBucket >   FunctionBucket_ptr;
  class SolverBucket;                                                // predeclaration
  typedef std::shared_ptr< SolverBucket >     SolverBucket_ptr;
  class FunctionalBucket;                                            // predeclaration
  typedef std::shared_ptr< FunctionalBucket > FunctionalBucket_ptr;
  class GenericDetectors;                                            // predeclaration
  typedef std::shared_ptr< GenericDetectors > GenericDetectors_ptr;
  class ReferencePoint;                                              // predeclaration
  typedef std::shared_ptr< ReferencePoint >   ReferencePoint_ptr;

  //*****************************************************************|************************************************************//
  // std shared pointers to dolfin objects
  //*****************************************************************|************************************************************//

  typedef std::shared_ptr< dolfin::GenericFunction >              GenericFunction_ptr;
  typedef std::shared_ptr< dolfin::Constant >                     Constant_ptr;
  typedef std::shared_ptr< dolfin::Expression >                   Expression_ptr;
  typedef std::shared_ptr< dolfin::Mesh >                         Mesh_ptr;
  typedef std::shared_ptr< const dolfin::Mesh >                   const_Mesh_ptr;
  typedef std::shared_ptr< dolfin::MeshFunction< std::size_t > >  MeshFunction_size_t_ptr;
  typedef std::shared_ptr< dolfin::MeshValueCollection< std::size_t > >  MeshValueCollection_size_t_ptr;
  typedef std::shared_ptr< dolfin::FunctionSpace >                FunctionSpace_ptr;
  typedef std::shared_ptr< dolfin::Function >                     Function_ptr;
  typedef std::shared_ptr< const dolfin::Function >               const_Function_ptr;
  typedef std::shared_ptr< dolfin::DirichletBC >                  DirichletBC_ptr;
  typedef std::shared_ptr< dolfin::Form >                         Form_ptr;
  typedef std::shared_ptr< dolfin::PETScMatrix >                  PETScMatrix_ptr;
  typedef std::shared_ptr< dolfin::PETScVector >                  PETScVector_ptr;
  typedef std::shared_ptr< const dolfin::PETScVector >            const_PETScVector_ptr;
  typedef std::shared_ptr< dolfin::GenericVector >                GenericVector_ptr;
  typedef std::shared_ptr< dolfin::XDMFFile >                     XDMFFile_ptr;
  typedef std::shared_ptr< dolfin::Array<double> >                Array_double_ptr;
  typedef std::shared_ptr< dolfin::SubDomain >                    SubDomain_ptr;
  typedef std::shared_ptr< dolfin::FunctionAssigner >             FunctionAssigner_ptr;

  //*****************************************************************|************************************************************//
  // iterators to std shared pointers in map pointer structures
  //*****************************************************************|************************************************************//

  typedef std::map< std::string, FunctionSpace_ptr >::iterator           FunctionSpace_it;
  typedef std::map< std::string, FunctionSpace_ptr >::const_iterator     FunctionSpace_const_it;
  typedef std::map< Mesh_ptr, FunctionSpace_ptr >::iterator              Mesh_FunctionSpace_it;
  typedef std::map< Mesh_ptr, FunctionSpace_ptr >::const_iterator        Mesh_FunctionSpace_const_it;
  typedef std::map< std::string, GenericFunction_ptr >::iterator         GenericFunction_it;
  typedef std::map< std::string, GenericFunction_ptr >::const_iterator   GenericFunction_const_it;
  typedef std::map< std::string, Expression_ptr >::iterator              Expression_it;
  typedef std::map< std::string, Expression_ptr >::const_iterator        Expression_const_it;
  typedef std::map< std::size_t, Expression_ptr >::iterator              size_t_Expression_it;
  typedef std::map< std::size_t, Expression_ptr >::const_iterator        size_t_Expression_const_it;
  typedef std::map< Mesh_ptr, XDMFFile_ptr >::iterator                   MeshXDMF_it;
  typedef std::map< Mesh_ptr, XDMFFile_ptr >::const_iterator             MeshXDMF_const_it;
  typedef std::map< std::string, bool_ptr >::iterator                    bool_ptr_it;
  typedef std::map< std::string, bool_ptr >::const_iterator              bool_ptr_const_it;

  struct om_key_seq{};
  struct om_key_hash{};

  template<typename Tkey, typename Tvalue>
  struct om_item {
    Tkey   first; 
    Tvalue second;
    om_item (Tkey &key, Tvalue &value) : first(key), second(value) {}
  };

  template<typename Tkey, typename Tvalue>
  using ordered_map = boost::multi_index::multi_index_container<
                            om_item<Tkey, Tvalue>,
                            boost::multi_index::indexed_by<
                              boost::multi_index::hashed_unique<boost::multi_index::tag<om_key_hash>, 
                                                                boost::multi_index::member<om_item<Tkey,Tvalue>, Tkey, &om_item<Tkey,Tvalue>::first>>,
                              boost::multi_index::sequenced<boost::multi_index::tag<om_key_seq> >
                            >
                          >;

  typedef boost::multi_index::index<ordered_map<const std::string,SystemBucket_ptr>,om_key_hash>::type::iterator       SystemBucket_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SystemBucket_ptr>,om_key_hash>::type::const_iterator SystemBucket_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SystemBucket_ptr>,om_key_seq>::type::iterator        SystemBucket_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SystemBucket_ptr>,om_key_seq>::type::const_iterator  SystemBucket_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,SolverBucket_ptr>,om_key_hash>::type::iterator       SolverBucket_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SolverBucket_ptr>,om_key_hash>::type::const_iterator SolverBucket_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SolverBucket_ptr>,om_key_seq>::type::iterator        SolverBucket_it;
  typedef boost::multi_index::index<ordered_map<const std::string,SolverBucket_ptr>,om_key_seq>::type::const_iterator  SolverBucket_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,FunctionBucket_ptr>,om_key_hash>::type::iterator       FunctionBucket_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionBucket_ptr>,om_key_hash>::type::const_iterator FunctionBucket_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionBucket_ptr>,om_key_seq>::type::iterator        FunctionBucket_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionBucket_ptr>,om_key_seq>::type::const_iterator  FunctionBucket_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,FunctionalBucket_ptr>,om_key_hash>::type::iterator       FunctionalBucket_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionalBucket_ptr>,om_key_hash>::type::const_iterator FunctionalBucket_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionalBucket_ptr>,om_key_seq>::type::iterator        FunctionalBucket_it;
  typedef boost::multi_index::index<ordered_map<const std::string,FunctionalBucket_ptr>,om_key_seq>::type::const_iterator  FunctionalBucket_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,Mesh_ptr>,om_key_hash>::type::iterator       Mesh_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Mesh_ptr>,om_key_hash>::type::const_iterator Mesh_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Mesh_ptr>,om_key_seq>::type::iterator        Mesh_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Mesh_ptr>,om_key_seq>::type::const_iterator  Mesh_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,MeshFunction_size_t_ptr>,om_key_hash>::type::iterator       MeshFunction_size_t_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,MeshFunction_size_t_ptr>,om_key_hash>::type::const_iterator MeshFunction_size_t_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,MeshFunction_size_t_ptr>,om_key_seq>::type::iterator        MeshFunction_size_t_it;
  typedef boost::multi_index::index<ordered_map<const std::string,MeshFunction_size_t_ptr>,om_key_seq>::type::const_iterator  MeshFunction_size_t_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,std::string>,om_key_hash>::type::iterator       string_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::string>,om_key_hash>::type::const_iterator string_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::string>,om_key_seq>::type::iterator        string_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::string>,om_key_seq>::type::const_iterator  string_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,GenericDetectors_ptr>,om_key_hash>::type::iterator       GenericDetectors_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,GenericDetectors_ptr>,om_key_hash>::type::const_iterator GenericDetectors_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,GenericDetectors_ptr>,om_key_seq>::type::iterator        GenericDetectors_it;
  typedef boost::multi_index::index<ordered_map<const std::string,GenericDetectors_ptr>,om_key_seq>::type::const_iterator  GenericDetectors_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,Form_ptr>,om_key_hash>::type::iterator       Form_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Form_ptr>,om_key_hash>::type::const_iterator Form_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Form_ptr>,om_key_seq>::type::iterator        Form_it;
  typedef boost::multi_index::index<ordered_map<const std::string,Form_ptr>,om_key_seq>::type::const_iterator  Form_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,DirichletBC_ptr>,om_key_hash>::type::iterator       DirichletBC_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,DirichletBC_ptr>,om_key_hash>::type::const_iterator DirichletBC_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,DirichletBC_ptr>,om_key_seq>::type::iterator        DirichletBC_it;
  typedef boost::multi_index::index<ordered_map<const std::string,DirichletBC_ptr>,om_key_seq>::type::const_iterator  DirichletBC_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,ReferencePoint_ptr>,om_key_hash>::type::iterator       ReferencePoint_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,ReferencePoint_ptr>,om_key_hash>::type::const_iterator ReferencePoint_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,ReferencePoint_ptr>,om_key_seq>::type::iterator        ReferencePoint_it;
  typedef boost::multi_index::index<ordered_map<const std::string,ReferencePoint_ptr>,om_key_seq>::type::const_iterator  ReferencePoint_const_it;

  typedef boost::multi_index::index<ordered_map<const std::string,std::vector<std::string>>,om_key_hash>::type::iterator       vector_string_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::vector<std::string>>,om_key_hash>::type::const_iterator vector_string_const_hash_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::vector<std::string>>,om_key_seq>::type::iterator        vector_string_it;
  typedef boost::multi_index::index<ordered_map<const std::string,std::vector<std::string>>,om_key_seq>::type::const_iterator  vector_string_const_it;

}

#endif
