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

import buckettools.systembucket
import buckettools.spud
import libspud
import os

class SpudSystemBucket(buckettools.systembucket.SystemBucket):
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space) 
     plus all the information necessary to populate that class using libspud."""
  
  def fill(self, optionpath, bucket):
    """Fill a system class with data describing that system using libspud and the given optionpath."""
    self.name       = libspud.get_option(optionpath+"/name")
    self.optionpath = optionpath
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol").split(os.linesep)[0]
    self.bucket     = bucket

    self.mesh_name        = libspud.get_option(optionpath+"/mesh/name")
    mesh_optionpath  = "/geometry/mesh::"+self.mesh_name
    self.cell             = libspud.get_option(mesh_optionpath+"/source/cell")

    self.fields = []
    for j in range(libspud.option_count(optionpath+"/field")):
      field_optionpath = optionpath+"/field["+`j`+"]"
      field = buckettools.spud.SpudFunctionBucket()
      # get all the information about this field from the options dictionary
      field.fill(field_optionpath, self, j)
      # let the system know about this field
      self.fields.append(field)
      # remove the local copy of this field
      del field

    self.coeffs = []
    for j in range(libspud.option_count(optionpath+"/coefficient")):
      coeff_optionpath = optionpath+"/coefficient["+`j`+"]"
      coeff = buckettools.spud.SpudFunctionBucket()
      # get all the information about this coefficient from the options dictionary
      coeff.fill(coeff_optionpath, self, j)
      # let the system know about this coefficient
      self.coeffs.append(coeff)
      # remove the local copy of this coefficient
      del coeff
    
    self.special_coeffs = []
    if libspud.have_option("/timestepping"):
      coeff_optionpath = "/timestepping/timestep/coefficient::Timestep"
      coeff = buckettools.spud.SpudFunctionBucket()
      # get all the information about this coefficient from the options dictionary
      coeff.fill(coeff_optionpath, self, 0)
      # let the system know about this coefficient
      self.special_coeffs.append(coeff)
      # remove the local copy of this coefficient
      del coeff
    
    self.solvers = []
    for j in range(libspud.option_count(optionpath+"/nonlinear_solver")):
      solver_optionpath = optionpath+"/nonlinear_solver["+`j`+"]"
      solver = buckettools.spud.SpudSolverBucket()
      # get all the information about this nonlinear solver from the options dictionary
      solver.fill(solver_optionpath, self)
      # let the system know about this solver
      self.solvers.append(solver)
      # done with this nonlinear solver
      del solver

    self.functionals = []
    for j in range(libspud.option_count(optionpath+"/functional")):
      functional_optionpath = optionpath+"/functional["+`j`+"]"
      functional = buckettools.spud.SpudFunctionalBucket()
      # get all the information about this functional from the options dictionary
      functional.fill(functional_optionpath, self)
      # let the system know about this functional
      self.functionals.append(functional)
      # done with this functional
      del functional


