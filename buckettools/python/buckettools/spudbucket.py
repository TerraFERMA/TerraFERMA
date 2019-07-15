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
import buckettools.bucket
import buckettools.spud

class SpudBucket(buckettools.bucket.Bucket):
  def fill(self):
    """Fill a bucket class with data describing a set of mixedfunctionspace systems using libspud, the given optionpath."""

    self.systems = []

    parameters_optionpath = "/global_parameters/ufl"
    if libspud.have_option(parameters_optionpath):
      self.parameters = libspud.get_option(parameters_optionpath)

    # loop over the systems in the options tree
    for i in range(libspud.option_count("/system")):
      system_optionpath = "/system["+repr(i)+"]"
      system = buckettools.spud.SpudSystemBucket()
      # get all the information about this system from the options dictionary
      system.fill(system_optionpath, self)
      # let the bucket know about this system
      self.systems.append(system)
      # done with this system
      del system

