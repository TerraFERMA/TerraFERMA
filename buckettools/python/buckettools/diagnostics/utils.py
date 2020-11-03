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
Some useful utility functions
"""

import copy
import time
import unittest

import buckettools.diagnostics.optimise as optimise

def CanLen(input):
  """
  Return whether it is valid to call len on the supplied input
  """
  
  try:
    len(input)
    return True
  except TypeError:
    return False
    
def ExpandList(input):
  """
  Return a list equal to the input with all sub-elements expanded
  """

  parentType = input.__class__
  outputList = []
  for item in input:
    # Trap the case that comparing the item with the input produces an array (as
    # is the case for numpy arrays)
    try:
      equalsInput = item == input
      if CanLen(equalsInput):
        if len(equalsInput) == 1:
          equalsInput = equalsInput[0]
        else:
          equalsInput = False
    except ValueError:
      equalsInput = False
      
    if (not CanLen(equalsInput) or len(equalsInput) == 1) and equalsInput:
      # Trap the case that iterating through the input produces the input (as
      # is the case for single character strings)
      outputList.append(item)
    # If the item is of the same type as the input expand it, unless the item
    # is a string and the input is not
    elif isinstance(item, parentType) or (CanLen(item) and not isinstance(item, str)):
      for subitem in ExpandList(item):
        outputList.append(subitem)
    else:
      outputList.append(item)
  
  return outputList
    
def FormLine(inputList, delimiter = " ", newline = True):
  """
  Form a delimited line out of the supplied list, (by default) terminated with a
  newline
  """
  
  inputList = ExpandList(inputList)

  line = ""
  for i in range(len(inputList)):
    line += str(inputList[i])
    if i < len(inputList) - 1:
      line += delimiter
  if newline:
    line += "\n"
  
  return line
  
# Jumping through Python 2.3 hoops for cx1 - in >= 2.4, can just pass a cmp
# argument to list.sort
class Sorter:
  def __init__(self, key, value):
    self._key = key
    self._value = value
    
    return
    
  def GetKey(self):
    return self._key
  
  def GetValue(self):
    return self._value
    
  def __cmp__(self, val):
    if self._key > val:
      return 1
    elif self._key == val:
      return 0
    else:
      return -1
  
def KeyedSort(keys, *values, **kargs):
  """
  Return a re-ordering of values, with the remapping equal to the sort mapping
  of keys. Each key must correspond to a unique value. Optional keyword
  arguments:
    returnSortedKeys  Return the sorted keys as well as the sorted values
  """
  
  returnSortedKeys = False
  for key in kargs:
    if key == "returnSortedKeys":
      returnSortedKeys = True
    else:
      raise Exception("Unexpected keyword argument \"" + key + "\"")
  
  sortList = []
  for i, key in enumerate(keys):
    sortList.append(Sorter(key, tuple([subValues[i] for subValues in values])))
    
  sortList.sort()
  if optimise.DebuggingEnabled():
    for i in range(len(sortList) - 1):
      if sortList[i].GetKey() == sortList[i + 1].GetKey():
        assert(sortList[i].GetValue() == sortList[i].GetValue())
  
  result = []
  if returnSortedKeys:
    result.append([sorter.GetKey() for sorter in sortList])
  
  for i in range(len(values)):
    result.append([sorter.GetValue()[i] for sorter in sortList])
  
  if len(result) == 1:
    result = result[0]
  else:
    result = tuple(result)
  
  return result
  
