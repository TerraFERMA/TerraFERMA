#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

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


import argparse
try:
  import argcomplete
except ImportError:
  pass
from copy import copy
import dolfin
import sys
import subprocess
import numpy

parser = argparse.ArgumentParser( \
                       description="""This takes a dolfin .xml mesh file """ +\
                       """and extracts the submesh with the given region id(s).""")
parser.add_argument('filename', metavar='filename', type=str,
                    help='specify the name of the dolfin .xml or .xml.gz file')
parser.add_argument('-r', '--regionids', action='store', metavar='split_ids', dest='regionids', type=int, nargs='+', required=True,
                    help='specify the region ids where the submesh is to be extracted')
parser.add_argument('-o', '--outputfilename', action='store', metavar='filename', dest='outputfilename', type=str, default=None, required=False,
                    help='specify the output filename (defaults to the input filename_submesh)')
try:
  argcomplete.autocomplete(parser)
except NameError:
  pass
args = parser.parse_args()

# check that the filename ends with the right format
if args.filename[-4:]!=".xml" and args.filename[-7:]!=".xml.gz":
    sys.stderr.write("Mesh filename must end in .xml or .xml.gz.\n")
    parser.print_help()
    sys.exit(1)

# organize the filename
fullname = args.filename
if fullname[-4:]==".xml":
  basename = fullname[:-4]
  extension = ".xml"
else:
  basename = fullname[:-7]
  extension = ".xml.gz"

outname = args.outputfilename
if outname is None: outname = basename+"_submesh"
if outname[-4:]==".xml":
  outname = outname[:-4]
  extension = ".xml"
elif outname[-7:]==".xml.gz":
  outname = outname[:-7]
  extension = ".xml.gz"

# read in the mesh
mesh = dolfin.Mesh(fullname)
print args.regionids
submesh = dolfin.SubMesh(mesh, args.regionids)

submesh_file_out = dolfin.File(outname+".xml")
submesh_file_out << submesh
if extension[-3:]==".gz":
  subprocess.call(['gzip', '-f', outname+".xml"])


