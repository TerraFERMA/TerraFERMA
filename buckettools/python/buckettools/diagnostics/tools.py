# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# (originally part of fluidity, modified for buckettools by Cian Wilson)

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

"""
TerraFERMA related tools
"""

import glob
import os
import subprocess
import tempfile
import unittest
import numpy

import buckettools.statfile as statfile

import buckettools.diagnostics.debug as debug
import buckettools.diagnostics.calc as calc
import buckettools.diagnostics.utils as utils
  
class Stat:
  """
  Class for handling stat files. Similiar to the dictionary returned by
  parser, but with some annoying features fixed.
  """

  def __init__(self, filename = None, delimiter = "%", includeMc = False, subsample = 1):
    self.SetDelimiter(delimiter)
    self._s = {}
    if not filename is None:
      self.Read(filename, includeMc = includeMc, subsample = subsample)
      
    return
    
  def __getitem__(self, key):
    """
    Index into the .stat with the given key (or path)
    """
  
    def SItem(s, key, delimiter):
      if key in s:
        return s[key]
      else:
        keySplit = key.split(delimiter)
        for i in range(len(keySplit)):
          key = utils.FormLine(keySplit[:i], delimiter = delimiter, newline = False)
          if key in s:
            if isinstance(s[key], dict):
              try:
                # Tolerate a failure when recursing, as the key may have been
                # eroneously split
                return SItem(s[key], utils.FormLine(keySplit[i:], delimiter = delimiter, newline = False), delimiter)
              except Exception:
                pass
            else:
              return s[key]

        raise Exception("Key not found")
    
    item = SItem(self._s, key, self._delimiter)
    if isinstance(item, dict):
      subS = Stat(delimiter = self._delimiter)
      subS._s = item
      return subS
    else:
      return item
    
  def __setitem__(self, key, value):
    keySplit = self.SplitPath(key)
    assert(len(keySplit) > 0)
    s = self
    for key in keySplit[:-1]:
      if s.haskey(key):
        s = s[key]
        assert(isinstance(s, Stat))
      else:
        s._s[key] = {}
        s = s[key]
    s._s[keySplit[-1]] = value
    
    return
    
  def __str__(self):
    paths = self.Paths()
    string = "Stat file:\n"
    for i, path in enumerate(paths):
      string += path
      if i < len(paths):
        string += "\n"
        
    return string
    
  def _PathSplit(self, path):
    """
    Return a list of keys into the stat dictionary for the supplied path
    """
    
    def SPathSplit(s, delimiter, path):
      pathSplit = path.split(delimiter)
      index = 0
      newPath = pathSplit[index]
      while not newPath in s.keys():
        index += 1
        newPath += delimiter + pathSplit[index]
        
      paths = []
      paths.append(newPath)
      if isinstance(s[newPath], dict):
        paths += SPathSplit(s[newPath], delimiter, path[len(newPath) + len(delimiter):])
        
      return paths
  
    return SPathSplit(self._s, self._delimiter, path)
    
  def haskey(self, key):
    return self.HasPath(key)
    
  def keys(self):
    return self.Paths()
    
  def GetDelimiter(self):
    return self._delimiter 
   
  def SetDelimiter(self, delimiter):
    self._delimiter = delimiter
    
    return
    
  def Paths(self):
    """
    Return all valid paths
    """
  
    def SPaths(s, delimiter, base = ""):
      if len(base) > 0:
        base += delimiter
    
      paths = []
      for key in s.keys():
        if isinstance(s[key], dict):
          paths += SPaths(s[key], delimiter, base = base + key)
        else:
          paths.append(base + key)
          
      return paths
    
    return SPaths(self._s, self._delimiter)
    
  def PathLists(self):
    """
    Return all valid paths as a series of key lists
    """
    
    def SPathLists(s, delimiter, base = []):    
      paths = []
      for key in s.keys():
        if isinstance(s[key], dict):
          paths += SPathLists(s[key], delimiter, base = base + [key])
        else:
          paths.append(base + [key])
          
      return paths
    
    return SPathLists(self._s, self._delimiter)
    
  def HasPath(self, path):
    """
    Return whether the supplied path is valid for this Stat
    """
    
    try:
      self[path]
      return True
    except Exception:
      return False
    
  def FormPath(self, *args):
    path = ""
    for i, arg in enumerate(args):
      path += arg
      if i < len(args) - 1:
        path += self._delimiter
    
    return path
    
  def FormPathFromList(self, pathList):
    path = ""
    for i, entry in enumerate(pathList):
      path += entry
      if i < len(pathList) - 1:
        path += self._delimiter
        
    return path
    
  def SplitPath(self, path):
    return path.split(self._delimiter)
    
  def Read(self, filename, includeMc = False, subsample = 1):
    """
    Read a .stat file
    """
    
    def ParseRawS(s, delimiter):    
      newS = {}
      for key1 in s.keys():
        assert(not key1 in ["val", "value"])
        if isinstance(s[key1], dict):
          if len(s[key1].keys()) == 1 and s[key1].keys()[0] in ["val", "value"]:
            newS[str(key1)] = s[key1][s[key1].keys()[0]]
          else:
            subS = ParseRawS(s[key1], delimiter)
            newS[str(key1)] = {}
            for key2 in subS.keys():
              newS[str(key1)][str(key2)] = subS[key2]
        else:        
          rank = len(s[key1].shape)
          if rank > 1:
            assert(rank == 2)
            if includeMc:
              # Add in this vector
              
              # parser gives this in an inconvenient matrix order. Take the
              # transpose here to make life easier.
              newS[str(key1)] = s[key1].transpose()
              
            # Add in the vector field components
            for i in range(len(s[key1])):
              newS[str(key1) + delimiter + str(i + 1)] = s[key1][i]
          else:
            try:
              # Add in this scalar
              newS[str(key1)] = s[key1]
            except TypeError:
              debug.deprint("Type error for data " + str(s[key1]), 0)
              raise Exception("ParseRawS failure")
            except ValueError:
              debug.deprint("Value error for data " + str(s[key1]), 0)
              raise Exception("ParseRawS failure")
          
      return newS
        
    debug.dprint("Reading .stat file: " + filename)
    if subsample == 1:
      # Handle this case separately, as it's convenient to be backwards
      # compatible
      statParser = statfile.parser(filename)
    else:
      statParser = statfile.parser(filename, subsample = subsample)

    self._s = ParseRawS(statParser, self._delimiter)
    
    if "ElapsedTime" in self.keys():
      t = self["ElapsedTime"]
      if t.shape[0] > 0:
        debug.dprint("Time range: " + str((t[0], t[-1])))
      else:
        debug.dprint("Time range: No data")
    
    return
    
def JoinStat(*args):
  """
  Joins a series of stat files together. Useful for combining checkpoint .stat
  files. Selects data in later stat files over earlier stat files. Assumes
  data in stat files are sorted by ElapsedTime.
  """

  nStat = len(args)
  assert(nStat > 0)
  times = [stat["ElapsedTime"] for stat in args]
  
  startT = [t[0] for t in times]
  permutation = utils.KeyedSort(startT, range(nStat))
  stats = [args[index] for index in permutation]
  startT = [startT[index] for index in permutation]
  times = [times[index] for index in permutation]
  
  endIndices = numpy.array([len(time) for time in times], dtype = int)
  for i, t in enumerate(times[:-1]):
    for j, time in enumerate(t):
      if calc.AlmostEquals(startT[i + 1], time, tolerance = 1.0e-6):
        endIndices[i] = max(j - 1, 0)
        break
      elif startT[i + 1] < time:
        endIndices[i] = j
        break
  debug.dprint("Time ranges:")
  if len(times) > 0:
    for i in range(nStat): 
      debug.dprint((startT[i], times[i][endIndices[i] - 1]))
  else:
    debug.dprint("No data")
    
  dataIndices = numpy.empty(len(args) + 1, dtype = int)
  dataIndices[0] = 0
  for i, index in enumerate(endIndices):
    dataIndices[i + 1] = dataIndices[i] + index
    
  stat = stats[0]
  data = {}
  for key in stat.keys():
    arr = stat[key]
    shape = list(arr.shape)
    shape[0] = dataIndices[-1]
    data[key] = numpy.empty(shape, dtype = arr.dtype)
    data[key][:dataIndices[1]] = arr[:endIndices[0]]
    data[key][dataIndices[1]:] = calc.Nan()
  delimiter = stat.GetDelimiter()
  
  for i in range(1, nStat):
    stat = stats[i]
    for key in stat.keys(): 
      arr = stat[key]
      if not key in data:
        shape = list(arr.shape)
        shape[0] = dataIndices[-1]
        data[key] = numpy.empty(shape, dtype = arr.dtype)
        data[key][:dataIndices[i]] = calc.Nan()
        data[key][dataIndices[i + 1]:] = calc.Nan()
      data[key][dataIndices[i]:dataIndices[i + 1]] = arr[:endIndices[i]]
  
  output = Stat(delimiter = delimiter)
  for key in data.keys():
    output[key] = numpy.array(data[key])
  
  return output
                
