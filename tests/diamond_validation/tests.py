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


#from optparse import OptionParser
import sys
import glob
import os
import xml.dom.minidom
import xml.parsers.expat

import diamond.debug as debug
import diamond.schema as schema

debug.SetDebugLevel(0)

#optparser=OptionParser(usage='usage: %prog <basedirectoryname>',
#                       add_help_option=True,
#                       description="""This script tests options files validity against schemas.""")

#(options, argv) = optparser.parse_args()

#if len(argv)==0:
#    basedirectoryname = "."
#else:
#    basedirectoryname = argv[0]

ignoredOptionsFiles = []
try:
  warningsFile = open("ignored", "r")
except IOError:
  warningsFile = None
if warningsFile is None:
  ignoredOptionsFiles = []
else:
  ignoredOptionsFiles = warningsFile.readlines()
  for i, filename in enumerate(ignoredOptionsFiles):
    ignoredOptionsFiles[i] = filename[:-1]

class DiamondValidator:
  def __init__(self, rootDir):
    self._rootDir = rootDir
    self.Reset()
    
    return
    
  def Reset(self):
    self._passes = 0
    self._optionErrors = {}
    
    return
    
  def ValidateOptionsFiles(self, schemafile, testDir, depth, extension = None, xmlRootNode = None):
    debug.dprint("Validating options file against schema: " + schemafile, 0)
  
    schemafile = os.path.join(self._rootDir, schemafile)
    sch = schema.Schema(schemafile)

    if not extension is None:
      debug.dprint("Testing files with extension: " + extension, 0)
      for filename in self._TestFiles(extension, testDir, depth): 
        optionsTree = sch.read(filename)
        lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
        if len(lost_eles) + len(added_eles) + len(lost_attrs) + len(added_attrs) == 0 and optionsTree.valid:
          debug.dprint(filename + " : Pass", 0)
          self._passes += 1
        else:
          debug.dprint(filename + " : Fail", 0)
          self._optionErrors[filename] = (lost_eles, added_eles, lost_attrs, added_attrs)
          
    if not xmlRootNode is None:
      debug.dprint("Testing xml files with root node: " + xmlRootNode, 0)
      for filename in self._TestFiles("xml", testDir, depth):
        try:
          xmlParse = xml.dom.minidom.parse(filename)
        except xml.parsers.expat.ExpatError:
          continue
        rootEles = xmlParse.getElementsByTagName(xmlRootNode)
        if len(rootEles) == 0:
          continue
        optionsTree = sch.read(filename)
        lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
        if len(lost_eles) + len(added_eles) + len(lost_attrs) + len(added_attrs) == 0 and optionsTree.valid:
          debug.dprint(filename + " : Pass", 0)
          self._passes += 1
        else:
          debug.dprint(filename + " : Fail", 0)
          self._optionErrors[filename] = (lost_eles, added_eles, lost_attrs, added_attrs)
    
    return
    
  def _TestFiles(self, extension, testDir, depth):
    filenames = []
    baseDir = os.path.join(self._rootDir, testDir)
    for i in range(depth + 1):
      filenames += glob.glob(os.path.join(baseDir, "*." + extension))
      baseDir = os.path.join(baseDir, "*")
    
    return filenames
    
  def Passes(self):
    return self._passes
    
  def OptionErrors(self):
    return self._optionErrors

validator = DiamondValidator(rootDir = os.path.join(os.path.pardir, os.path.pardir))

validator.ValidateOptionsFiles(schemafile = os.path.join("buckettools", "schemas", "terraferma.rng"), testDir = "tests", depth = 1, extension = "tfml")

passes = validator.Passes()
optionErrors = validator.OptionErrors()
    
print "Summary of options files with failures:"
for filename in optionErrors.keys():
  print filename
print "Passes: " + str(passes)
print "Failures: " + str(len(optionErrors))

def test_validity_failures():
  failures = []
  for filename in optionErrors.keys():
    if not filename in ignoredOptionsFiles:
      failures.append(filename)
  print "Summary of options files with failures:"
  for filename in failures:
    print filename
  print "Failures: " + str(len(failures))
  assert(len(failures) == 0)

