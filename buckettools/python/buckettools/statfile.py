# (originally part of fluidity, modified for buckettools by Cian Wilson)

# Copyright (C) 2006 Imperial College London and others.
# 
# Please see the AUTHORS file in the main source directory for a full list
# of copyright holders.
#
# Prof. C Pain
# Applied Modelling and Computation Group
# Department of Earth Science and Engineering
# Imperial College London
#
# amcgsoftware@imperial.ac.uk
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation,
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

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

import array
import os
import re
import numpy
from xml.dom.minidom import parseString
from xml.dom.minidom import Document

class creator(dict):
  """Class to create stat files. The stat entries are defined using  
   creator[system][name][statistic] 
   or
   creator[name][statistic].

   Constants can be added with the add_constant function.

   Example:
     from buckettools.statfile import creator
     c=creator("my_stat.stat")
     c.add_constant({"time": 1.0})
     c[('Material1', 'Speed', 'max')] = 1.0
     c.write()
  """

  def __init__(self, filename):
    self.filename = filename
    self.initialised = False
    self.constants = {}

  def add_constant(self, constant):
    if self.initialised:
      print("Constant can only be added before the the first write() call")
      return
    self.constants.update(constant)

  def write(self):
    if not self.initialised:
      f = open(self.filename, "w")
      # Create the minidom document
      doc = Document()
      # Create the <header> element
      header = doc.createElement("header")
      doc.appendChild(header)
      self.header = [] # We save the header for verification before every write_stat
      # Write the constants
      for const_k, const_v in list(self.constants.items()):
        const_element = doc.createElement("constant")
        const_element.setAttribute("name", str(const_k))
        const_element.setAttribute("type", "string")
        const_element.setAttribute("value", str(const_v))
        header.appendChild(const_element)
      # Create the stat elements
      column = 1
      for stat_k, stat_v in list(self.items()):
        stat_element = doc.createElement("field")
        stat_element.setAttribute("column", str(column))
        if len(stat_k) == 2:
          stat_element.setAttribute("name", stat_k[0])
          stat_element.setAttribute("statistic", stat_k[1])
        elif len(stat_k) == 3:
          stat_element.setAttribute("system", stat_k[0])
          stat_element.setAttribute("name", stat_k[1])
          stat_element.setAttribute("statistic", stat_k[2])
        else:
          print("Element ", stat_k, " must have length 2 or 3")
          exit()
        if isinstance(stat_v, (list, numpy.ndarray)):
          stat_element.setAttribute("components", str(len(stat_v)))
          column += len(stat_v)
        else:
          column += 1
        header.appendChild(stat_element)
        self.header.append(stat_k)
      self.initialised = True
      try:
            f.write(doc.toprettyxml(indent="  "))
      finally:
            f.close()
      # Now call the write function again to actually write the first values
      self.write() 
      return
    # Here the header is written and we only want to append data. So lets load the file in append mode 
    f = open (self.filename, "a")
    # Check that the dictionary and the header are consistent 
    if set(self) != set(self.header):
      print("Error: Columns may not change after initialisation of the stat file.")
      print("Columns you are trying to write: ", self)
      print("Columns in the header: ", self.header)
      exit()
    output = ""
    for stat_k in self.header:
      stat_v = self[stat_k]
      if isinstance(stat_v, (list, numpy.ndarray)):
        for v in stat_v:
          output = output + "  " + str(v)
      else:
        output = output + "  " + str(stat_v)
    output = output+os.linesep  
    try:
          f.write(output)
    finally:
          f.close()

class parser(dict):
    """Parse a .stat file. The resulting mapping object is a hierarchy
of dictionaries. Most entries are of the form:

   parser[system][field][statistic].

for example:

   p=stat_parser(filename)
   p['Material1']['Speed']['max']
"""

    def __init__(self, filename, subsample = 1):
    
      assert(subsample > 0)

      statfile=open(filename, "r")
      header_re=re.compile(r"</header>")
      xml="" # xml header.

      # extract the xml header stopping when </header> is reached.
      while 1:
        line=statfile.readline()
        if line=="":
          raise Exception("Unable to read .stat file header")
        xml=xml+line
        if re.search(header_re, line):
          break

      # now parse the xml.
      parsed=parseString(xml)
      
      binaryFormat = False
      constantEles = parsed.getElementsByTagName("constant")
      self.constants = {}
      for ele in constantEles:
        name = ele.getAttribute("name")
        vtype = ele.getAttribute("type")
        value = ele.getAttribute("value")
        self.constants[name] = value
        if name == "format":
          assert(vtype == "string")
          if value == "binary":
            binaryFormat = True
       
      nColumns = 0  
      for field in parsed.getElementsByTagName("field"):
        components = field.getAttribute("components")
        if components:
          nColumns += int(components)
        else:
          nColumns += 1

      if binaryFormat:           
        for ele in constantEles:
          name = ele.getAttribute("name")
          vtype = ele.getAttribute("type")
          value = ele.getAttribute("value")
          if name == "real_size":
            assert(vtype == "integer")            
            real_size = int(value)            
            if real_size == 4:
              realFormat = 'f'
            elif real_size == 8:
              realFormat = 'd'
            else:
              raise Exception("Unexpected real size: " + str(real_size))
          elif name == "integer_size":
            assert(vtype == "integer")            
            integer_size = int(value)            
            if not integer_size == 4:
              raise Exception("Unexpected integer size: " + str(integer_size))
       
        nOutput = int( (os.path.getsize(filename + ".dat") / (nColumns * real_size)) / subsample)
        
        columns = numpy.empty((nColumns, nOutput))
        statDatFile = open(filename + ".dat", "rb")   
        index = 0     
        while True:
          values = array.array(realFormat)
          try:
            values.fromfile(statDatFile, nColumns)
          except EOFError:
            break
            
          for i, value in enumerate(values):
            columns[i][index] = value
            
          index += 1
          if index >= nOutput:
            # Ignore incomplete lines
            break
          if subsample > 1:
            # Ignore non-sampled lines
            statDatFile.seek(real_size * (subsample - 1) * nColumns, 1)        
        statDatFile.close()
        assert(index == nOutput)
      else:
        columns = [[] for i in range(nColumns)]
        lineNo = 0
        for line in statfile:
          entries = list(map(float, line.split()))
          # Ignore non-sampled lines
          if len(entries) == len(columns) and (lineNo % subsample) == 0:
            list(map(list.append, columns, entries))
          elif len(entries) != len(columns):
            raise Exception("Incomplete line %d: expected %d, but got %d columns" % (lineNo, len(columns), len(entries)))
          lineNo = lineNo + 1
        columns = numpy.array(columns)
              
      for field in parsed.getElementsByTagName("field"):
        system=field.getAttribute("system")
        name=field.getAttribute("name")
        column=field.getAttribute("column")
        statistic=field.getAttribute("statistic")
        components=field.getAttribute("components")

        if system:
          if system not in self:
            self[system]={}
          current_dict=self[system]
        else:
          current_dict=self

        if name not in current_dict:
          current_dict[name]={}

        if components:
            column=int(column)
            components=int(components)
            current_dict[name][statistic]=columns[column-1:column-1+components]
        else:
            current_dict[name][statistic]=columns[int(column)-1]

