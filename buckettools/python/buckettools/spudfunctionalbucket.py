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
import buckettools.functionalbucket
import os

class SpudFunctionalBucket(buckettools.functionalbucket.FunctionalBucket):
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl) 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, system, function=None):
    """Fill a functional class with data describing that functional using libspud, the given optionpath and the system it's based on."""
    try:
      self.name   = libspud.get_option(optionpath+"/name")
    except libspud.SpudKeyError:
      if function is None:
        self.name = ""
      else:
        self.name = function.name
    self.symbol   = libspud.get_option(optionpath+"/ufl_symbol").split(os.linesep)[0]
    self.form     = libspud.get_option(optionpath)+os.linesep
    self.system   = system
    self.function = function
    self.form_representation = libspud.get_option(optionpath+"/form_representation/name")
    if libspud.have_option(optionpath+"/quadrature_degree"):
      self.quadrature_degree = libspud.get_option(optionpath+"/quadrature_degree")
    self.quadrature_rule = libspud.get_option(optionpath+"/quadrature_rule/name")


  
