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
    self.checkpoint = False
    if self.type == "TimeSeries":
      self.grids = grid.findall("Grid")
      self.grids.sort(key=lambda g: float(g.find("Time").attrib["Value"]))
      self.times = np.asarray([float(g.find("Time").attrib["Value"]) for g in self.grids])
    elif self.type == "Grid":
      self.grids = [grid]
      self.times = np.asarray([0.0])
    else:
      # FIXME: just assume that the other format is a checkpoint file
      self.type = "Checkpoint"  # the type will be a field name so let's make up a type for this
      self.checkpoint = True
      self.grids = self.tree.findall("/Domain/Grid")
      # Assuming that the times are the same for all grids
      self.times = np.asarray([float(g.find("Time").attrib["Value"]) for g in self.grids[0].getchildren()])      

  def _getgrid(self, name=None, tindex=-1, time=None):
    """
    Return the grid at the given tindex or time. 
    If time is given it overrides tindex, which will be ignored.
    name is only used if this is a checkpoint xdmf, if it is not given then the first field grid is returned.
    """
    times = self.gettimes(name=name)
    if time is not None: tindex = np.abs(times-time).argmin()
    if tindex >= len(times) or (tindex < 0 and abs(tindex) > len(times)):
      raise IndexError("tindex out of range")

    if self.checkpoint:
      i = 0
      if name is not None:
        names = self.getfieldnames()
        try:
          i = names.index(name)
        except ValueError:
          self._unknownfieldname()
      fieldgrid = self.grids[i]
      grid = fieldgrid.getchildren()[tindex]
    else:
      grid = self.grids[tindex]
    return grid

  def _getattr(self, name, tindex=-1, time=None, index=-1):
    """
    Return the named attribute at the given tindex, time and index.
    If time is given it overrides tindex, which will be ignored.
    """
    grid = self._getgrid(name=name, tindex=tindex, time=time)
    attrs = grid.findall("Attribute[@Name='"+name+"']")
    # checkpoints should only have a single attribute so sorting isn't necessary
    if not self.checkpoint: attrs.sort(key=lambda a: int(a.find("DataItem").text.split("/")[-1]))
    if index >= len(attrs) or (index < 0 and abs(index) > len(attrs)):
      raise IndexError("attribute missing or index out of range")
    return attrs[index]

  def _getvalues(self, element):
    """
    Given an element from an xdmf file that contains a DataItem,
    return the numpy array of the data.
    """
    dataitems = element.findall("DataItem")
    if len(dataitems) > 1:
      # if there's more than one data item, take the float (only expect this to apply to checkpoints)
      dataitem = [di for di in dataitems if di.attrib.get('NumberType', 'Float')=='Float'][0]
    else:
      dataitem = dataitems[0]
    i = dataitem.text.index(":")
    h5filename = dataitem.text[:i]
    h5keys = [k for k in dataitem.text[i+1:].split("/") if k != '']
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

  def _unknownfieldname(self):
    """
    Raise an exception for an unknown field name.
    """
    raise Exception("Unknown field name.")

  def _unknowncelltype(self):
    """
    Raise an exception for an unknown VTK cell type.
    """
    raise Exception("Unknown VTK cell type.")

  def _unknownrank(self):
    """
    Raise an exception for an unknown function rank.
    """
    raise Exception("Unknown function rank.")

  def _untestedcelltype(self):
    """
    Warn about an untested VTK cell type.
    """
    warnings.warn("Untested VTK cell type.  Please test.", stacklevel=2)

  def getlocations(self, name=None, tindex=-1, time=None):
    """
    Return the coordinates of the grid point locations at the given tindex or time.

    If time is given it overrides tindex.
    
    'name' is only used if this is a checkpoint xdmf, if it is not given then the locations of the first field are returned.
    """
    grid = self._getgrid(name=name, tindex=tindex, time=time)
    if grid.find("Geometry") is None:
      grid = self._getgrid(name=name, tindex=0)
    geom = grid.find("Geometry")
    coords = self._getvalues(geom)
    return coords

  def gettopology(self, name=None, tindex=-1, time=None):
    """
    Return the topology of the grid points at the given tindex or time.
    
    If time is given it overrides tindex.

    'name' is only used if this is a checkpoint xdmf, if it is not given then the topology of the first field is returned.
    """
    
    grid = self._getgrid(name=name, tindex=tindex, time=time)
    if grid.find("Topology") is None:
      grid = self._getgrid(name=name, tindex=0)
    topo = grid.find("Topology")
    topoarr = self._getvalues(topo)
    return topoarr

  def gettimes(self, name=None):
    """
    Return the times available in the xdmf.

    'name' is only used if this is a checkpoint xdmf, if not set the first set of times is returned.
    """
    if self.checkpoint:
      i = 0
      if name is not None:
        names = self.getfieldnames()
        try:
          i = names.index(name)
        except ValueError:
          self._unknownfieldname()
      fieldgrid = self.grids[i]
      times = np.asarray([float(g.find("Time").attrib["Value"]) for g in fieldgrid.getchildren()])
    else:
      times = self.times
    return times

  def getnindices(self, name=None, tindex=-1, time=None):
    """
    Return the number of fields available at the given tindex or time.
    
    If time is given it overrides tindex.

    If name is not given it defaults to the first field name available.  All fields are assumed to have the same number of available
versions.

    If no fields (attributes) are found, 1 is returned.
    """
    grid = self._getgrid(name=name, tindex=tindex, time=time)
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
    return self._getvalues(attr)

  def getfieldnames(self, tindex=-1, time=None):
    """
    Return the field anmes at the given tindex or time.

    If time is given it overrides tindex.
    """
    if self.checkpoint:
      names = [g.attrib["Name"] for g in self.grids]
    else:
      grid = self._getgrid(tindex=tindex, time=time)
      attrs = grid.findall("Attribute")
      names = list(set([a.attrib["Name"] for a in attrs]))
    return names

  #def probedata(self, coordinates, name, tindex=-1, time=None, index=-1):
  #  """
  #  Interpolate field values at the given coordinates
  #  """

  def _vtucheckpoint(self, names=[], tindex=-1, time=None, index=-1):
    from buckettools import vtktools as vtk
    import dolfin as df
    from collections import OrderedDict
    vtu = vtk.vtu()

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

    functions = {}
    mdeg = 1
    disc = False
    cell = None
    for name in names:
      attr = self._getattr(name, tindex=tindex, time=time, index=index)
      family = attr.attrib["ElementFamily"]
      degree = int(attr.attrib["ElementDegree"])
      fcell  = attr.attrib["ElementCell"]
      frank  = attr.attrib["AttributeType"]
      if cell is None: cell = fcell
      if fcell != cell: raise Exception("Mixed cell types.")
      functions[name] = (family, degree, fcell, frank)
      mdeg = max(degree, mdeg)
      disc = disc or (family != "CG" and degree > 0)

    # make a mesh
    # FIXME: assuming the same mesh for all fields (so defaulting to first)
    coords = self.getlocations(tindex=tindex, time=time)
    topo = self.gettopology(tindex=tindex, time=time)
    d = coords.shape[-1]
    if cell=="interval": d = 1
    mesh = df.Mesh()
    editor = df.MeshEditor()
    editor.open(mesh, cell, d, d)
    editor.init_vertices(coords.shape[0])
    editor.init_cells(topo.shape[0])
    for i,c in enumerate(coords): editor.add_vertex(i, c[:d])
    for i,c in enumerate(topo): editor.add_cell(i, c)
    editor.close()

    family = "DG" if disc else "CG"
    degree = 2 if mdeg > 1 else 1
    sV = df.FunctionSpace(mesh, family, degree)
    vV = df.VectorFunctionSpace(mesh, family, degree)
    tV = df.TensorFunctionSpace(mesh, family, degree)
    sfunc = df.Function(sV)
    vfunc = df.Function(vV)
    tfunc = df.Function(tV)

    coordmap = OrderedDict()
    dofmap = sV.dofmap()
    element = sV.element()
    for cell in df.cells(mesh):
      dofcoord = element.tabulate_dof_coordinates(cell)
      celldofs = dofmap.cell_dofs(cell.index())
      for i,cd in enumerate(celldofs):
        coordmap[cd] = dofcoord[i]
    n = len(celldofs) # FIXME: assuming single cell type

    # some dofmaps
    # to file
    ldofmap = np.zeros(len(coordmap), dtype=int)
    for i,ci in enumerate(coordmap.keys()): ldofmap[ci] = i
    # from file
    rdofmap = np.zeros(len(coordmap), dtype=int)
    for i,ci in enumerate(ldofmap): rdofmap[ci] = i
    # cell-based
    cdofmap = np.asarray(list(range(topo.shape[0])), dtype=int)

    # add points
    points = vtk.vtk.vtkPoints()
    points.SetDataTypeToDouble()
    for c in coordmap.values():
      cp = c.tolist()+[0.0]*(3-d)
      points.InsertNextPoint(cp[0], cp[1], cp[2])
    vtu.ugrid.SetPoints(points)

    # get cell details
    # FIXME: assuming simplicial
    celltypes = {(0,1) : vtk.vtk.VTK_VERTEX,
                 (1,1) : vtk.vtk.VTK_LINE,
                 (1,2) : vtk.vtk.VTK_QUADRATIC_EDGE,
                 (2,1) : vtk.vtk.VTK_TRIANGLE,
                 (2,2) : vtk.vtk.VTK_QUADRATIC_TRIANGLE,
                 (3,1) : vtk.vtk.VTK_TETRA,
                 (3,2) : vtk.vtk.VTK_QUADRATIC_TETRA,
                }
    celltype = celltypes.get((d,degree)) # FIXME: assuming single cell type
    cellorder = list(range(n))
    if celltype is None:
      self._unknowncelltype()
    elif celltype == vtk.vtk.VTK_QUADRATIC_TRIANGLE:
      cellorder = [0,1,2,5,3,4]
    elif celltype == vtk.vtk.VTK_QUADRATIC_TETRA:
      cellorder = [0,1,2,3,9,6,8,7,5,4]
    
    dfxdmf = df.XDMFFile(self.filename)

    # add the cells
    for cell in df.cells(mesh):
      celldofs = dofmap.cell_dofs(cell.index())
      idList = vtk.vtk.vtkIdList()
      for i in range(n): idList.InsertNextId(ldofmap[celldofs[cellorder[i]]])
      vtu.ugrid.InsertNextCell(celltype, idList)

    # add the fields
    pointdata = vtu.ugrid.GetPointData()
    celldata = vtu.ugrid.GetCellData()
    for name in names:
      lfamily, ldegree, lcell, lrank = functions[name]

      if lrank == "Scalar":
        lfunc = sfunc
      elif lrank == "Vector":
        lfunc = vfunc
      elif lrank == "Tensor":
        lfunc = tfunc
      else:
        self._unknownrank()

      if lfamily == family and ldegree == degree:
        cfunc = lfunc
      else:
        if (lrank == "Scalar") or (lrank == "Vector" and lfamily in ["RT", "DRT", "BDM", "N1curl", "N2curl"]):
          V = df.FunctionSpace(mesh, lfamily, ldegree)
        elif lrank == "Vector":
          V = df.VectorFunctionSpace(mesh, lfamily, ldegree)
        else:
          V = df.TensorFunctionSpace(mesh, lfamily, ldegree)
        cfunc = df.Function(V)

      dfxdmf.read_checkpoint(cfunc, name, tindex)

      if (lfamily == "DG" and ldegree == 0) or (lfamily == family and ldegree == degree):
        lfunc = cfunc
      else:
        lfunc.interpolate(cfunc)

      dpp = 1
      fdim = np.asarray(lfunc.value_shape(), dtype=int).prod()
      indices = list(range(fdim))
      if lrank == "Vector":
        dpp = 3
      elif lrank == "Tensor":
        dpp = 9
        indices[2] = 3
        indices[3] = 4

      lrdofmap = rdofmap
      if lfamily == "DG" and ldegree == 0: lrdofmap = cdofmap
      ndofs = len(lrdofmap)

      vals = np.zeros((ndofs, dpp))
      lfuncs = [lfunc]
      if lrank != "Scalar": lfuncs = lfunc.split(deepcopy=True)
      for i in range(len(lfuncs)):
        vals[:, indices[i]] = lvec = lfuncs[i].vector()[lrdofmap]

      data = vtk.vtk.vtkDoubleArray()
      data.SetNumberOfValues(vals.size)
      if lrank != "Scalar": data.SetNumberOfComponents(dpp)
      data.SetName(name)
      vals = vals.flatten()
      for i, v in enumerate(vals): data.SetValue(i, v)
      if lfamily == "DG" and ldegree == 0:
        celldata.AddArray(data)
        celldata.SetActiveScalars(name)
      else:
        pointdata.AddArray(data)
        pointdata.SetActiveScalars(name)

    dfxdmf.close()

    return vtu

  def _vtutsgrid(self, names=[], tindex=-1, time=None, index=-1):
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
      vals = self._getvalues(attr)
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

  def vtu(self, names=[], tindex=-1, time=None, index=-1):
    """
    Convert the given tindex or time and index to a vtu.

    If time is given it overrides tindex.

    If names is given, only those fields are output.  If not supplied all fields are included.
    """
    vtu = None
    if self.checkpoint:
      vtu = self._vtucheckpoint(names=names,tindex=tindex,time=time,index=index)
    else:
      vtu = self._vtutsgrid(names=names,tindex=tindex,time=time,index=index)

    return vtu

