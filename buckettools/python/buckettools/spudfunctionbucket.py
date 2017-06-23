# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

import sys
import libspud
import buckettools.functionbucket
import os

class SpudFunctionBucket(buckettools.functionbucket.FunctionBucket):
  """A class that stores all the information necessary to write the ufl for a function (field or coefficient) 
     plus all the information necessary to populate that class using libspud.
     Note that the class has limited ufl production because much of this is system dependent."""

  def fill(self, optionpath, system, index):
    """Fill a function class with data describing that function using libspud, the given optionpath and the system its based on."""
    self.name       = libspud.get_option(optionpath+"/name")
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol").split(os.linesep)[0]
    self.system     = system
    self.index      = index
    self.type       = libspud.get_option(optionpath+"/type/name")

    self.rank   = libspud.get_option(optionpath+"/type/rank/name")
    self.family   = None
    self.degree = None
    self.functional = None
    if self.type != "Constant":
      self.family = libspud.get_option(optionpath+"/type/rank/element/family")
      self.degree = libspud.get_option(optionpath+"/type/rank/element/degree")

    self.size     = None
    self.shape    = None
    self.symmetry = None
    if self.rank == "Vector":
      if libspud.have_option(optionpath+"/type/rank/element/size"):
        self.size = libspud.get_option(optionpath+"/type/rank/element/size")
    elif self.rank == "Tensor":
      if libspud.have_option(optionpath+"/type/rank/element/shape"):
        self.shape = libspud.get_option(optionpath+"/type/rank/element/shape")
      if libspud.have_option(optionpath+"/type/rank/element/symmetric"):
        self.symmetry = True

    self.enrichment_family = None
    self.enrichment_degree = None
    if libspud.have_option(optionpath+"/type/rank/element/enrichment"):
      self.enrichment_family = libspud.get_option(optionpath+"/type/rank/element/enrichment/element/family")
      self.enrichment_degree = libspud.get_option(optionpath+"/type/rank/element/enrichment/element/degree")

    if libspud.have_option(optionpath+"/type/rank/element/quadrature_rule"):
      self.quadrature_rule = libspud.get_option(optionpath+"/type/rank/element/quadrature_rule/name")

    # this should be restricted by the schema to Constant coefficients:
    if libspud.have_option(optionpath+"/type/rank/value/functional"):
      functional_optionpath = optionpath+"/type/rank/value/functional"
      functional = buckettools.spud.SpudFunctionalBucket()
      functional.fill(functional_optionpath, self.system, self)
      self.functional = functional
    
    self.cpp = []
    for k in range(libspud.option_count(optionpath+"/type/rank/initial_condition")):
      cpp_optionpath = optionpath+"/type/rank/initial_condition["+`k`+"]"
      if libspud.have_option(cpp_optionpath+"/cpp"):
        cpp_name = libspud.get_option(cpp_optionpath+"/name")
        cppexpression = buckettools.spud.SpudCppExpressionBucket()
        # get all the information about this expression from the options dictionary
        cppexpression.fill(cpp_optionpath, cpp_name, self)
        # let the field know about this cpp expression
        self.cpp.append(cppexpression)
        # done with this expression
        del cppexpression

    for j in range(libspud.option_count(optionpath+"/type/rank/boundary_condition")):
      bc_optionpath = optionpath+"/type/rank/boundary_condition["+`j`+"]"
      bc_name = libspud.get_option(bc_optionpath+"/name")
      for k in range(libspud.option_count(bc_optionpath+"/sub_components")):
        bc_comp_optionpath = bc_optionpath+"/sub_components["+`k`+"]"
        bc_comp_name = libspud.get_option(bc_comp_optionpath+"/name")
        
        cpp_optionpath = bc_comp_optionpath+"/type"
        if libspud.have_option(cpp_optionpath+"/cpp"):
          cpp_name = bc_name + bc_comp_name
          
          cppexpression = buckettools.spud.SpudCppExpressionBucket()
          # get all the information about this expression from the options dictionary
          cppexpression.fill(cpp_optionpath, cpp_name, self)
          # let the field know about this cpp expression
          self.cpp.append(cppexpression)
          # done with this expression
          del cppexpression

    for k in range(libspud.option_count(optionpath+"/type/rank/value")):
      cpp_optionpath = optionpath+"/type/rank/value["+`k`+"]"
      if libspud.have_option(cpp_optionpath+"/cpp"):
        cpp_name = libspud.get_option(cpp_optionpath+"/name")
        cppexpression = buckettools.spud.SpudCppExpressionBucket()
        # get all the information about this expression from the options dictionary
        cppexpression.fill(cpp_optionpath, cpp_name, self)
        # let the field know about this cpp expression
        self.cpp.append(cppexpression)
        # done with this expression
        del cppexpression

