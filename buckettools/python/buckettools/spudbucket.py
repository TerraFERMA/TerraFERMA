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
import re
import os

class SpudBucket(buckettools.bucket.Bucket):
  def fill(self):
    """Fill a bucket class with data describing a set of mixedfunctionspace systems using libspud, the given optionpath."""

    self.meshes = {}
    # loop over the meshes in the options tree
    for i in range(libspud.option_count("/geometry/mesh")):
      mesh_optionpath = "/geometry/mesh["+repr(i)+"]"
      mesh_name = libspud.get_option(mesh_optionpath+"/name")
      self.meshes[mesh_name] = libspud.get_option(mesh_optionpath+"/source/cell")

    visualization_optionpath = "/io/visualization/element"
    try:
      self.viselementfamily = libspud.get_option(visualization_optionpath+"/family")
      self.viselementdegree = libspud.get_option(visualization_optionpath+"/degree")
    except libspud.SpudKeyError:
      self.viselementfamily = "CG"
      self.viselementdegree = 1

    parameters_optionpath = "/global_parameters/ufl"
    if libspud.have_option(parameters_optionpath):
      self.parameters = libspud.get_option(parameters_optionpath)

    self.systems = []
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

    self.cpplibraries = []
    cpplibraries_optionpath = "/global_parameters/cpp/libraries"
    if libspud.have_option(cpplibraries_optionpath):
      # regular expression to split up values string into a list
      # the list may be comma(,), semicolon(;), space ( ) or newline (\n) delimited
      #NODE                     EXPLANATION
      #--------------------------------------------------------------------------------
      #  (?:                      group, but do not capture (1 or more times
      #                           (matching the most amount possible)):
      #--------------------------------------------------------------------------------
      #    [^,; \n]                 any character except: ',', ';', ' ',
      #                             '\n' (newline)
      #--------------------------------------------------------------------------------
      #  )+                       end of grouping
      # - from http://rick.measham.id.au/paste/explain.pl
      r = re.compile(r'(?:[^,; '+os.linesep+'])+')
      self.cpplibraries = r.findall(libspud.get_option(cpplibraries_optionpath))

