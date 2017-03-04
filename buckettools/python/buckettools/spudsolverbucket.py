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

import libspud
import buckettools.solverbucket

class SpudSolverBucket(buckettools.solverbucket.SolverBucket):
  """A class that stores all the information necessary to write the ufl for a system of forms (i.e. linear or bilinear) associated with a solver 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, system):
    """Fill a solver class with data describing a system of forms using libspud, the given optionpath and the system its based on."""
    self.name = libspud.get_option(optionpath+"/name")
    newoptionpath = optionpath+"/type"
    self.type = libspud.get_option(newoptionpath+"/name")
     
    if libspud.have_option(newoptionpath+"/preamble"):
      self.preamble = libspud.get_option(newoptionpath+"/preamble")+"\n"

    self.form_names = []
    self.forms = []
    self.form_symbols = []
    self.form_ranks = []
    self.fill_subforms(newoptionpath)
    prefix = system.name+"_"+self.name+"_"
    self.fill_solverforms(newoptionpath, prefix=prefix)
    
    self.form_representation = libspud.get_option(newoptionpath+"/form_representation/name")

    if libspud.have_option(newoptionpath+"/quadrature_degree"):
      self.quadrature_degree = libspud.get_option(newoptionpath+"/quadrature_degree")

    self.quadrature_rule = libspud.get_option(newoptionpath+"/quadrature_rule/name")

    self.system = system

  def fill_subforms(self, optionpath, prefix=""):
    for i in range(libspud.option_count(optionpath+"/form")):
      form_optionpath = optionpath+"/form["+`i`+"]"
      self.form_names.append(prefix+libspud.get_option(form_optionpath+"/name"))
      self.forms.append(libspud.get_option(form_optionpath)+"\n")
      self.form_symbols.append(libspud.get_option(form_optionpath+"/ufl_symbol").split("\n")[0])
      self.form_ranks.append(int(libspud.get_option(form_optionpath+"/rank")))

  def fill_solverforms(self, optionpath, prefix=""):
    schurpc_optionpath = optionpath+"/composite_type::schur/schur_preconditioner::user"
    if (libspud.have_option(schurpc_optionpath)):
      self.fill_subforms(schurpc_optionpath, prefix=prefix)
    
    fs_optionpath = optionpath+"/fieldsplit"
    ls_optionpath = optionpath+"/linear_solver/preconditioner"
    pc_optionpath = optionpath+"/preconditioner"
    if (libspud.have_option(fs_optionpath)):
      for i in range(libspud.option_count(fs_optionpath)):
        newoptionpath = fs_optionpath+"["+`i`+"]"
        name = libspud.get_option(newoptionpath+"/name")
        self.fill_solverforms(newoptionpath+"/linear_solver/preconditioner", prefix=prefix+name+"_")
    elif (libspud.have_option(ls_optionpath)):
      self.fill_solverforms(ls_optionpath, prefix=prefix)
    elif (libspud.have_option(pc_optionpath)):
      self.fill_solverforms(pc_optionpath, prefix=prefix)


