import array
import exceptions
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
     from buckettools import stat_creator
     c=stat_creator("my_stat.stat")
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
      print "Constant can only be added before the the first write() call"
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
      for const_k, const_v in self.constants.items():
        const_element = doc.createElement("constant")
        const_element.setAttribute("name", str(const_k))
        const_element.setAttribute("type", "string")
        const_element.setAttribute("value", str(const_v))
        header.appendChild(const_element)
      # Create the stat elements
      column = 1
      for stat in self.keys():
        stat_element = doc.createElement("field")
        stat_element.setAttribute("column", str(column))
        if len(stat) == 2:
          stat_element.setAttribute("name", stat[0])
          stat_element.setAttribute("statistic", stat[1])
        elif len(stat) == 3:
          stat_element.setAttribute("system", stat[0])
          stat_element.setAttribute("name", stat[1])
          stat_element.setAttribute("statistic", stat[2])
        else:
          print "Element ", stat, " must have length 2 or 3"
          exit()
        header.appendChild(stat_element)
        self.header.append(stat)
        column = column+1
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
      print "Error: Columns may not change after initialisation of the stat file."
      print "Columns you are trying to write: ", self
      print "Columns in the header: ", self.header
      exit()
    output = ""
    for stat in self.header:
      output = output + "  " + str(self[stat])
    output = output + "\n"  
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

      statfile=file(filename, "r")
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
      for ele in constantEles:
        name = ele.getAttribute("name")
        type = ele.getAttribute("type")
        value = ele.getAttribute("value")
        if name == "format":
          assert(type == "string")
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
          type = ele.getAttribute("type")
          value = ele.getAttribute("value")
          if name == "real_size":
            assert(type == "integer")            
            real_size = int(value)            
            if real_size == 4:
              realFormat = 'f'
            elif real_size == 8:
              realFormat = 'd'
            else:
              raise Exception("Unexpected real size: " + str(real_size))
          elif name == "integer_size":
            assert(type == "integer")            
            integer_size = int(value)            
            if not integer_size == 4:
              raise Exception("Unexpected integer size: " + str(real_size))
       
        nOutput = (os.path.getsize(filename + ".dat") / (nColumns * real_size)) / subsample
        
        columns = numpy.empty((nColumns, nOutput))
        statDatFile = file(filename + ".dat", "rb")   
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
          entries = map(float, line.split())
          # Ignore non-sampled lines
          if len(entries) == len(columns) and (lineNo % subsample) == 0:
            map(list.append, columns, entries)
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
          if not self.has_key(system):
            self[system]={}
          current_dict=self[system]
        else:
          current_dict=self

        if not current_dict.has_key(name):
          current_dict[name]={}

        if components:
            column=int(column)
            components=int(components)
            current_dict[name][statistic]=columns[column-1:column-1+components]
        else:
            current_dict[name][statistic]=columns[int(column)-1]

