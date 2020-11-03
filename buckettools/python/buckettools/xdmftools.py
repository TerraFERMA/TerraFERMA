from lxml import etree
import h5py
import sys
import numpy as np
import os
import warnings

def h5find(h5, keys):
  """
  Recurse through a nested h5 file using a list of keys.
  """
  if len(keys) > 1:
    return h5find(h5[keys[0]], keys[1:])
  elif len(keys) == 1:
    return h5[keys[0]]
  else:
    raise Exception("Unknown h5 keys type")

class XDMF(object):
  def __init__(self, filename):
    """
    A XDMF class with the given filename.
    """
    self.filename = filename
    self.path = os.path.split(filename)[0]
    self.tree = etree.parse(filename)
    grid = self.tree.find("//Grid")
    self.type = grid.attrib["Name"]
    if self.type not in ["TimeSeries", "Grid"]: 
      self._unknownxdmftype()
    if self.type == "TimeSeries":
      self.grids = grid.findall("Grid")
      self.grids.sort(key=lambda g: float(g.find("Time").attrib["Value"]))
      self.times = np.asarray([float(g.find("Time").attrib["Value"]) for g in self.grids])
    elif self.type == "Grid":
      self.grids = [grid]
      self.times = np.asarray([0.0])

  def _getgrid(self, tindex=-1, time=None):
    """
    Return the grid at the given tindex or time. 
    If time is given it overrides tindex, which will be ignored.
    """
    if time is not None: tindex = np.abs(self.times-time).argmin()
    if tindex >= len(self.grids) or (tindex < 0 and abs(tindex) > len(self.grids)):
      raise IndexError("tindex out of range")
    return self.grids[tindex]

  def _getattr(self, name, tindex=-1, time=None, index=-1):
    """
    Return the named attribute at the given tindex, time and index.
    If time is given it overrides tindex, which will be ignored.
    """
    grid = self._getgrid(tindex=tindex, time=time)
    attrs = grid.findall("Attribute[@Name='"+name+"']")
    attrs.sort(key=lambda a: int(a.find("DataItem").text.split("/")[-1]))
    if index >= len(attrs) or (index < 0 and abs(index) > len(attrs)):
      raise IndexError("attribute missing or index out of range")
    return attrs[index]

  def _getarray(self, element):
    """
    Given an element from an xdmf file that contains a DataItem,
    return the numpy array of the data.
    """
    dataitem = element.find("DataItem")
    i = dataitem.text.index(":")
    h5filename = dataitem.text[:i]
    h5keys = dataitem.text[i+2:].split("/")
    h5file = h5py.File(os.path.join(self.path, h5filename), "r")
    dataset = h5find(h5file, h5keys)
    array = np.asarray(dataset)
    h5file.close()
    return array

  def _unknownxdmftype(self):
    """
    Raise an exception for an unknown XDMF file type.
    """
    raise Exception("Unknown XDMF type.")

  def _unknowncelltype(self):
    """
    Raise an exception for an unknown VTK cell type.
    """
    raise Exception("Unknown VTK cell type.")

  def _untestedcelltype(self):
    """
    Warn about an untested VTK cell type.
    """
    warnings.warn("Untested VTK cell type.  Please test.", stacklevel=2)

  def getlocations(self, tindex=-1, time=None):
    """
    Return the coordinates of the grid point locations at the given tindex or time.

    If time is given it overrides tindex.
    """
    if self.type in ["TimeSeries", "Grid"]:
      grid = self._getgrid(tindex=tindex, time=time)
      if grid.find("Geometry") is None:
        grid = self._getgrid(tindex=0)
      geom = grid.find("Geometry")
    else:
      self._unknownxdmftype()
    coords = self._getarray(geom)
    return coords

  def gettopology(self, tindex=-1, time=None):
    """
    Return the topology of the grid points at the given tindex or time.
    
    If time is given it overrides tindex.
    """
    if self.type in ["TimeSeries", "Grid"]:
      grid = self._getgrid(tindex=tindex, time=time)
      if grid.find("Topology") is None:
        grid = self._getgrid(tindex=0)
      topo = grid.find("Topology")
    else:
      self._unknownxdmftype()
    topoarr = self._getarray(topo)
    return topoarr

  def gettimes(self):
    """
    Return the times available in the xdmf.
    """
    return self.times

  def getnindices(self, name=None, tindex=-1, time=None):
    """
    Return the number of fields available at the given tindex or time.
    
    If time is given it overrides tindex.

    If name is not given it defaults to the first field name available.  All fields are assumed to have the same number of available
versions.

    If no fields (attributes) are found, 1 is returned.
    """
    grid = self._getgrid(tindex=tindex, time=time)
    if name is None:
      # assuming that nindices is the same for all fields in a given timestep
      names = self.getfieldnames(tindex=tindex, time=time)
      if len(names) > 0: name = names[0]
    if name is None:
      return 1
    else:
      attrs = grid.findall("Attribute[@Name='"+name+"']")
      return len(attrs)

  def getfield(self, name, tindex=-1, time=None, index=-1):
    """
    Return the named fields at the given tindex or time and index.

    If time is given it overrides tindex.
    """
    attr = self._getattr(name, tindex=tindex, time=time, index=index)
    return self._getarray(attr)

  def getfieldnames(self, tindex=-1, time=None):
    """
    Return the field anmes at the given tindex or time.

    If time is given it overrides tindex.
    """
    grid = self._getgrid(tindex=tindex, time=time)
    attrs = grid.findall("Attribute")
    return list(set([a.attrib["Name"] for a in attrs]))

  #def probedata(self, coordinates, name, tindex=-1, time=None, index=-1):
  #  """
  #  Interpolate field values at the given coordinates
  #  """

  def vtu(self, names=[], tindex=-1, time=None, index=-1):
    """
    Convert the given tindex or time and index to a vtu.

    If time is given it overrides tindex.

    If names is given, only those fields are output.  If not supplied all fields are included.
    """
    if self.type not in ['TimeSeries', 'Grid']:
      self._unknownxdmftype()

    from buckettools import vtktools as vtk
    
    vtu = vtk.vtu()

    # add points
    points = vtk.vtk.vtkPoints()
    points.SetDataTypeToDouble()
    coords = self.getlocations(tindex=tindex, time=time)
    d = coords.shape[-1]
    for c in coords:
      cp = c.tolist()+[0.0]*(3-d)
      points.InsertNextPoint(cp[0], cp[1], cp[2])
    vtu.ugrid.SetPoints(points)

    # get cell details
    topo = self.gettopology(tindex=tindex, time=time)
    n = topo.shape[-1]
    grid = self._getgrid(tindex=tindex, time=time)
    if grid.find("Topology") is None:
      grid = self._getgrid(tindex=0)
    topotype = grid.find("Topology").attrib.get("TopologyType")
    celltypes = {"PolyVertex"    : vtk.vtk.VTK_VERTEX,
                 "PolyLine"      : vtk.vtk.VTK_LINE,
                 "Edge_3"        : vtk.vtk.VTK_QUADRATIC_EDGE,
                 "Triangle"      : vtk.vtk.VTK_TRIANGLE,
                 "Triangle_6"    : vtk.vtk.VTK_QUADRATIC_TRIANGLE,
                 "Quadrilateral" : vtk.vtk.VTK_QUAD,
                 "Quad_8"        : vtk.vtk.VTK_QUADRATIC_QUAD,
                 "Tetrahedron"   : vtk.vtk.VTK_TETRA,
                 "Tet_10"        : vtk.vtk.VTK_QUADRATIC_TETRA,
                 "Hexahedron"    : vtk.vtk.VTK_HEXAHEDRON,
                 "Hex_20"        : vtk.vtk.VTK_QUADRATIC_HEXAHEDRON
                }
    celltype = celltypes.get(topotype) # assuming single cell type
    cellorder = list(range(n))
    if celltype is None:
      self._unknowncelltype()
    elif celltype == vtk.vtk.VTK_QUADRATIC_TRIANGLE:
      cellorder = [0,1,2,5,3,4]
    elif celltype == vtk.vtk.VTK_QUADRATIC_TETRA:
      cellorder = [0,1,2,3,9,6,8,7,5,4]
    
    # add the cells
    for t in topo:
      idList = vtk.vtk.vtkIdList()
      for i in cellorder:
        idList.InsertNextId(t[i])
      vtu.ugrid.InsertNextCell(celltype, idList)

    # add the fields
    if names == []:
      names = self.getfieldnames(tindex=tindex, time=time)
    else:
      if isinstance(names, str):
        names = [names]
      else:
        # check if iterable
        try:
          iter(names)
        except TypeError:
          raise Exception("field names must be a string or list of strings")
    pointdata = vtu.ugrid.GetPointData()
    celldata = vtu.ugrid.GetCellData()
    for name in names:
      attr = self._getattr(name, tindex=tindex, time=time, index=index)
      center = attr.attrib['Center']
      attrtype = attr.attrib['AttributeType']
      vals = self._getarray(attr)
      data = vtk.vtk.vtkDoubleArray()
      data.SetNumberOfValues(vals.size)
      if attrtype != 'Scalar': 
        data.SetNumberOfComponents(np.asarray(vals.shape[1:]).prod())
      data.SetName(name)
      vals = vals.flatten()
      for i, v in enumerate(vals):
        data.SetValue(i, v)
      if center == 'Node':
        pointdata.AddArray(data)
        pointdata.SetActiveScalars(name)
      elif center == 'Cell':
        celldata.AddArray(data)
        celldata.SetActiveScalars(name)
      else:
        raise Exception("Unknown data center.")

    return vtu

