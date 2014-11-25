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

import glob
import os
import sys
import argparse
import dolfin
from lxml import etree

parser = argparse.ArgumentParser( \
                       description="""A very simple script to plot DOLFIN meshfunctions.""")
parser.add_argument('file', action='store', metavar='file', type=str, nargs='+',
                    help='specify a filename or expression')
parser.add_argument('-m', '--mesh', metavar='mesh', action='store', type=str, dest='mesh', default=None, 
                    required=True,
                    help='specify the name of the mesh file')
parser.add_argument('-r', '--recursive', metavar='depth', action='store', type=int, dest='recurse', nargs='?', default=None, 
                    required=False, const=-1,
                    help='recursively search the directory tree for files (if no depth is specified full recursion will be used)')
args = parser.parse_args()

filenames = []
for f in args.file:
  
  if args.recurse is None:
    for filename in glob.glob(f):
      ext = filename.split('.')[-1]
      if ext == "xml":
        filenames.append(os.path.join(os.curdir, filename))
      else:
        print "Don't know how to deal with extension "+ext+".  Only know about xmls."
        sys.exit(1)
  else:
    
    if os.path.isabs(f):
      dirname = os.path.dirname(f)
    else:
      dirname = os.curdir
    dirname = os.path.normpath(dirname)
  
    for root, dirnames, files in os.walk(dirname, topdown=True):
      depth = string.count(root, os.path.sep)
      for filename in glob.glob1(os.path.join(root, os.path.dirname(f)), os.path.split(f)[-1]):
        ext = filename.split('.')[-1]
        if ext == "xml":
          filenames.append(os.path.join(os.path.join(root, os.path.dirname(f)), filename))
        else:
          print "Don't know how to deal with extension "+ext+".  Only know about xmls."
          sys.exit(1)
      if depth == args.recurse:
        dirnames[:] = []

mesh = dolfin.Mesh(args.mesh)

for filename in filenames:
  tree = etree.parse(filename)
  mvcs = tree.findall("//mesh_value_collection")
  if (len(mvcs)==0):
    print "No mesh_value_collection element found in %s, skipping."%(filename)
    continue
  elif (len(mvcs)>1):
    print "Multiple mesh_value_collection elements found in %s, skipping."%(filename)
    continue
  mftype = mvcs[0].attrib['type']
  mf = dolfin.MeshFunction(mftype, mesh, filename)
  mf.rename(filename, filename)
  dolfin.plot(mf)

dolfin.interactive()

