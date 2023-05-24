import math
import sys
import numpy
import vtk
import os

# All returned arrays are cast into either numpy or numarray arrays
arr=numpy.array

class vtu(object):
  """Unstructured grid object to deal with VTK unstructured grids."""
  def __init__(self, filename = None, tindex = -1, time = None, index = -1, family=None, degree=None):
    """Creates a vtu object by reading the specified file."""
    if filename is None:
      self.ugrid = vtk.vtkUnstructuredGrid()
    else:
      gridreader = None
      ext = os.path.splitext(filename)[-1]

      # "deal" with pvds
      if ext == ".pvd":
        from lxml import etree
        tree = etree.parse(filename)
        vtus = [element for element in tree.getroot().iterdescendants(tag="DataSet")]
        times = numpy.asarray([float(vtu.attrib["timestep"]) for vtu in vtus])
        filenames = [vtu.attrib["file"] for vtu in vtus]
        ltindex = tindex
        if time is not None: ltindex = numpy.abs(times-time).argmin()
        filename = filenames[ltindex]
        ext = os.path.splitext(filename)[-1]

      if ext == ".xdmf":
        from buckettools import xdmftools
        xdmf = xdmftools.XDMF(filename)
        self.ugrid = xdmf.vtu(tindex=tindex, time=time, index=index, family=family, degree=degree).ugrid
      else:
        if ext == ".vtu":
          gridreader=vtk.vtkXMLUnstructuredGridReader()
        elif ext == ".pvtu":
          gridreader=vtk.vtkXMLPUnstructuredGridReader()
        else:
          raise Exception("ERROR: don't recognise file extension" + filename)
        gridreader.SetFileName(filename)
        gridreader.Update()
        self.ugrid=gridreader.GetOutput()
      if self.ugrid.GetNumberOfPoints() + self.ugrid.GetNumberOfCells() == 0:
          raise Exception("ERROR: No points or cells found after loading vtu " + filename)
    self.filename=filename

  def GetScalarField(self, name):
    """Returns an array with the values of the specified scalar field."""
    try:
      pointdata=self.ugrid.GetPointData()
      vtkdata=pointdata.GetScalars(name)
      vtkdata.GetNumberOfTuples()
    except:
      try:
        celldata=self.ugrid.GetCellData()
        vtkdata=celldata.GetScalars(name)
        vtkdata.GetNumberOfTuples()
      except:
        raise Exception("ERROR: couldn't find point or cell scalar field data with name "+name+" in file "+self.filename+".")
    return arr([vtkdata.GetTuple1(i) for i in range(vtkdata.GetNumberOfTuples())])

  def GetScalarRange(self, name):
    """Returns the range (min, max) of the specified scalar field."""
    try:
      pointdata=self.ugrid.GetPointData()
      vtkdata=pointdata.GetScalars(name)
      vtkdata.GetRange()
    except:
      try:
        celldata=self.ugrid.GetCellData()
        vtkdata=celldata.GetScalars(name)
        vtkdata.GetRange()
      except:
        raise Exception("ERROR: couldn't find point or cell scalar field data with name "+name+" in file "+self.filename+".")
    return vtkdata.GetRange()

  def GetVectorField(self, name):
    """Returns an array with the values of the specified vector field."""
    try:
      pointdata=self.ugrid.GetPointData()
      vtkdata=pointdata.GetScalars(name)
      vtkdata.GetNumberOfTuples()
    except:
      try:
        celldata=self.ugrid.GetCellData()
        vtkdata=celldata.GetScalars(name)
        vtkdata.GetNumberOfTuples()
      except:
        raise Exception("ERROR: couldn't find point or cell vector field data with name "+name+" in file "+self.filename+".")
    return arr([vtkdata.GetTuple3(i) for i in range(vtkdata.GetNumberOfTuples())])

  def GetVectorNorm(self, name):
    """Return the field with the norm of the specified vector field."""
    v = self.GetVectorField(name)
    n = []

    try:
      from scipy.linalg import norm
    except ImportError:
      def norm(v):
        r = 0.0
        for x in v:
          r = r + x**2
        r = math.sqrt(r)
        return r

    for node in range(self.ugrid.GetNumberOfPoints()):
      n.append(norm(v[node]))

    return arr(n)

  def GetField(self,name):
    """Returns an array with the values of the specified field."""
    try:
      pointdata=self.ugrid.GetPointData()
      vtkdata=pointdata.GetArray(name)
      vtkdata.GetNumberOfTuples()
    except:
      try:
        celldata=self.ugrid.GetCellData()
        vtkdata=celldata.GetArray(name)
        vtkdata.GetNumberOfTuples()
      except:
        raise Exception("ERROR: couldn't find point or cell field data with name "+name+" in file "+self.filename+".")
    nc=vtkdata.GetNumberOfComponents()
    nt=vtkdata.GetNumberOfTuples()
    array=arr([vtkdata.GetValue(i) for i in range(nc*nt)])
    if nc==9:
      return array.reshape(nt,3,3)
    elif nc==4:
      return array.reshape(nt,2,2)
    else:
      return array.reshape(nt,nc)

  def GetFieldRank(self, name):
    """
    Returns the rank of the supplied field.
    """
    try:
      pointdata=self.ugrid.GetPointData()
      vtkdata=pointdata.GetArray(name)
      vtkdata.GetNumberOfTuples()
    except:
      try:
        celldata=self.ugrid.GetCellData()
        vtkdata=celldata.GetArray(name)
        vtkdata.GetNumberOfTuples()
      except:
        raise Exception("ERROR: couldn't find point or cell field data with name "+name+" in file "+self.filename+".")
    comps = vtkdata.GetNumberOfComponents()
    if comps == 1:
      return 0
    elif comps in [2, 3]:
      return 1
    elif comps in [4, 9]:
      return 2
    else:
      raise Exception("Field rank > 2 encountered")

  def Write(self, filename=None):
    """Writes the grid to a vtu file.

    If no filename is specified it will use the name of the file originally
    read in, thus overwriting it!
    """
    if filename is None:
      filename=self.filename
    if filename is None:
      raise Exception("No file supplied")
    if filename.endswith('pvtu'):
      gridwriter=vtk.vtkXMLPUnstructuredGridWriter()
    else:
      gridwriter=vtk.vtkXMLUnstructuredGridWriter()

    gridwriter.SetFileName(filename)
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      gridwriter.SetInput(self.ugrid)
    else:
      gridwriter.SetInputData(self.ugrid)
    gridwriter.Write()

  def AddScalarField(self, name, array):
    """Adds a scalar field with the specified name using the values from the array."""
    data = vtk.vtkDoubleArray()
    data.SetNumberOfValues(len(array))
    data.SetName(name)
    for i in range(len(array)):
      data.SetValue(i, array[i])

    if len(array) == self.ugrid.GetNumberOfPoints():
      pointdata=self.ugrid.GetPointData()
      pointdata.AddArray(data)
      pointdata.SetActiveScalars(name)
    elif len(array) == self.ugrid.GetNumberOfCells():
      celldata=self.ugrid.GetCellData()
      celldata.AddArray(data)
      celldata.SetActiveScalars(name)
    else:
      raise Exception("Length neither number of nodes nor number of cells")

  def AddVectorField(self, name, array):
    """Adds a vector field with the specified name using the values from the array."""
    n=array.size
    data = vtk.vtkDoubleArray()
    data.SetNumberOfComponents(array.shape[1])
    data.SetNumberOfValues(n)
    data.SetName(name)
    for i in range(n):
      data.SetValue(i, array.reshape(n)[i])

    if array.shape[0]==self.ugrid.GetNumberOfPoints():
      pointdata=self.ugrid.GetPointData()
      pointdata.AddArray(data)
      pointdata.SetActiveVectors(name)
    elif array.shape[0]==self.ugrid.GetNumberOfCells():
      celldata=self.ugrid.GetCellData()
      celldata.AddArray(data)
    else:
      raise Exception("Length neither number of nodes nor number of cells")

  def AddField(self, name, array):
    """Adds a field with arbitrary number of components under the specified name using."""
    n=array.size
    sh=arr(array.shape)
    data = vtk.vtkDoubleArray()
    # number of tuples is sh[0]
    # number of components is the product of the rest of sh
    data.SetNumberOfComponents(sh[1:].prod())
    data.SetNumberOfValues(n)
    data.SetName(name)
    flatarray=array.reshape(n)
    for i in range(n):
      data.SetValue(i, flatarray[i])

    if sh[0]==self.ugrid.GetNumberOfPoints():
      pointdata=self.ugrid.GetPointData()
      pointdata.AddArray(data)
    elif sh[0]==self.ugrid.GetNumberOfCells():
      celldata=self.ugrid.GetCellData()
      celldata.AddArray(data)
    else:
      raise Exception("Length neither number of nodes nor number of cells")

  def ApplyProjection(self, projection_x, projection_y, projection_z):
    """Applys a projection to the grid coordinates. This overwrites the existing values."""
    npoints = self.ugrid.GetNumberOfPoints ()
    for i in range (npoints):
      (x,y,z) = self.ugrid.GetPoint (i)
      new_x = eval (projection_x)
      new_y = eval (projection_y)
      new_z = eval (projection_z)
      self.ugrid.GetPoints ().SetPoint (i, new_x, new_y, new_z)

  def ApplyCoordinateTransformation(self,f):
    """Applys a coordinate transformation to the grid coordinates. This overwrites the existing values."""
    npoints = self.ugrid.GetNumberOfPoints ()
    
    for i in range (npoints):
      (x,y,z) = self.ugrid.GetPoint (i)
      newX = f(arr([x,y,z]),t=0)
      self.ugrid.GetPoints ().SetPoint (i, newX[0], newX[1], newX[2])

  def ApplyEarthProjection(self):
    """ Assume the input geometry is the Earth in Cartesian geometry and project to longatude, latitude, depth."""
    npoints = self.ugrid.GetNumberOfPoints ()
    
    earth_radius = 6378000.0
    rad_to_deg = 180.0/math.pi
    deg_to_rad = math.pi/180.0

    for i in range (npoints):
      (x,y,z) = self.ugrid.GetPoint (i)

      r = math.sqrt(x*x+y*y+z*z)
      depth = r - earth_radius
      longitude = rad_to_deg*math.atan2(y, x)
      latitude = 90.0 - rad_to_deg*math.acos(z/r)   
      
      self.ugrid.GetPoints ().SetPoint (i, longitude, latitude, depth)

  def ProbeData(self, coordinates, name):
    """Interpolate field values at these coordinates."""
    probe = VTUProbe(self.ugrid, coordinates)
    return probe.GetField(name)

  def RemoveField(self, name):
    """Removes said field from the unstructured grid."""
    pointdata=self.ugrid.GetPointData()
    pointdata.RemoveArray(name)

  def GetLocations(self):
    """Returns an array with the locations of the nodes."""
    vtkPoints = self.ugrid.GetPoints()
    if vtkPoints is None:
      vtkData = vtk.vtkDoubleArray()
    else:
      vtkData = vtkPoints.GetData()
    return arr([vtkData.GetTuple3(i) for i in range(vtkData.GetNumberOfTuples())])

  def GetCellPoints(self, id):
    """Returns an array with the node numbers of each cell (ndglno)."""
    idlist=vtk.vtkIdList()
    self.ugrid.GetCellPoints(id, idlist)
    return arr([idlist.GetId(i) for i in range(idlist.GetNumberOfIds())])

  def GetFieldNames(self):
    """Returns the names of the available fields."""
    vtkdata=self.ugrid.GetPointData()
    fieldnames = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    vtkdata=self.ugrid.GetCellData()
    fieldnames += [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    return fieldnames

  def GetPointCells(self, id):
    """Return an array with the elements which contain a node."""
    idlist=vtk.vtkIdList()
    self.ugrid.GetPointCells(id, idlist)
    return arr([idlist.GetId(i) for i in range(idlist.GetNumberOfIds())])

  def GetPointPoints(self, id):
    """Return the nodes connecting to a given node."""
    cells = self.GetPointCells(id)
    lst = []
    for cell in cells:
      lst = lst + list(self.GetCellPoints(cell))

    s = set(lst) # remove duplicates
    return arr(list(s)) # make into a list again

  def GetDistance(self, x, y):
    """Return the distance in physical space between x and y."""
    posx = self.ugrid.GetPoint(x)
    posy = self.ugrid.GetPoint(y)
    return math.sqrt(sum([(posx[i] - posy[i])**2 for i in range(len(posx))]))

  def Crop(self, min_x, max_x, min_y, max_y, min_z, max_z):
    """Trim off the edges defined by a bounding box."""
    trimmer = vtk.vtkExtractUnstructuredGrid()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      trimmer.SetInput(self.ugrid)
    else:
      trimmer.SetInputData(self.ugrid)
    trimmer.SetExtent(min_x, max_x, min_y, max_y, min_z, max_z)
    trimmer.Update()
    trimmed_ug = trimmer.GetOutput()

    self.ugrid = trimmed_ug

  def IntegrateField(self, field):
    """
    Integrate the supplied scalar field, assuming a linear representation on a
    tetrahedral mesh. Needs numpy-izing for speed.
    """

    assert field[0].shape in [(), (1,)]

    integral = 0.0

    n_cells = self.ugrid.GetNumberOfCells()
    vtkGhostLevels = self.ugrid.GetCellData().GetArray("vtkGhostLevels")
    for cell_no in range(n_cells):
      integrate_cell = True
      
      if vtkGhostLevels:
        integrate_cell = (vtkGhostLevels.GetTuple1(cell_no) == 0)

      if integrate_cell:
        Cell = self.ugrid.GetCell(cell_no)
        
        Cell_points = Cell.GetPoints ()
        nCell_points = Cell.GetNumberOfPoints()
        if nCell_points == 4:
          Volume = abs(Cell.ComputeVolume(Cell_points.GetPoint(0), \
                                            Cell_points.GetPoint(1), \
                                            Cell_points.GetPoint(2), \
                                            Cell_points.GetPoint(3)))
        elif nCell_points == 3:
          Volume = abs(Cell.TriangleArea(Cell_points.GetPoint(0), \
                                            Cell_points.GetPoint(1), \
                                            Cell_points.GetPoint(2)))
        else:
          raise Exception("Unexpected number of points: " + str(nCell_points))

        Cell_ids = Cell.GetPointIds()
        
        for point in range(Cell_ids.GetNumberOfIds()):
          PointId = Cell_ids.GetId(point)
          integral = integral + (Volume*field[PointId] / float(nCell_points))
          
    return integral

  def GetCellVolume(self, id):
    cell = self.ugrid.GetCell(id)
    pts = cell.GetPoints()
    if isinstance(cell, vtk.vtkTriangle):
      return cell.TriangleArea(pts.GetPoint(0), pts.GetPoint(1), pts.GetPoint(2))
    elif cell.GetNumberOfPoints() == 4:
      return abs(cell.ComputeVolume(pts.GetPoint(0), pts.GetPoint(1), pts.GetPoint(2), pts.GetPoint(3)))
    elif cell.GetNumberOfPoints() == 3:
      return abs(cell.ComputeVolume(pts.GetPoint(0), pts.GetPoint(1), pts.GetPoint(2)))
    else:
      raise Exception("Unexpected number of points")

  def GetFieldIntegral(self, name):
    """
    Integrate the named field.
    """

    return self.IntegrateField(self.GetField(name))

  def GetFieldRms(self, name):
    """
    Return the rms of the supplied scalar or vector field.
    """

    field = self.GetField(name)
    rank = self.GetFieldRank(name)
    if rank == 0:
      normField = arr([field[i] ** 2.0 for i in range(len(field))])
    elif rank == 1:
      normField = self.GetVectorNorm(name)
    else:
      raise Exception("Cannot calculate norm field for field rank > 1")
    volField = arr([1.0 for i in range(len(field))])
    rms = self.IntegrateField(normField)
    rms /= self.IntegrateField(volField)
    rms = numpy.sqrt(rms)

    return float(rms)

  def StructuredPointProbe(self, nx, ny, nz, bounding_box=None):
    """ Probe the unstructured grid dataset using a structured points dataset. """

    probe = vtk.vtkProbeFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      probe.SetSource(self.ugrid)
    else:
      probe.SetSourceData(self.ugrid)

    sgrid = vtk.vtkStructuredPoints()

    bbox = [0.0,0.0, 0.0,0.0, 0.0,0.0]
    if bounding_box==None:
      bbox = self.ugrid.GetBounds()
    else:
      bbox = bounding_box

    sgrid.SetOrigin([bbox[0], bbox[2], bbox[4]])

    sgrid.SetDimensions(nx, ny, nz)

    spacing = [0.0, 0.0, 0.0]
    if nx>1: spacing[0] = (bbox[1]-bbox[0])/(nx-1.0)
    if ny>1: spacing[1] = (bbox[3]-bbox[2])/(ny-1.0)
    if nz>1: spacing[2] = (bbox[5]-bbox[4])/(nz-1.0)

    sgrid.SetSpacing(spacing)

    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      probe.SetInput(sgrid)
    else:
      probe.SetInputData(sgrid)
    probe.Update ()

    return probe.GetOutput()

  ### Field manipulation methods ###

  def ManipulateField(self, fieldName, manipFunc, newFieldName = None):
    """
    Generic field manipulation method. Applies the supplied manipulation function
    manipFunc to the field fieldName. manipFunc must be a function of the form:

      def manipFunc(field, index):
        # ...
        return fieldValAtIndex
    """

    field = self.GetField(fieldName)
    if newFieldName is None or fieldName == newFieldName:
      self.RemoveField(fieldName)
      newFieldName = fieldName

    field = arr([manipFunc(field, i) for i in range(len(field))])
    self.AddField(newFieldName, field)

    return

  def AddFieldToField(self, fieldName, array, newFieldName = None):
    def ManipFunc(field, index):
      return field[index] + array[index]

    self.ManipulateField(fieldName, ManipFunc, newFieldName)

    return

  def SubFieldFromField(self, fieldName, array, newFieldName = None):
    def ManipFunc(field, index):
      return field[index] - array[index]

    self.ManipulateField(fieldName, ManipFunc, newFieldName)

    return

  def DotFieldWithField(self, fieldName, array, newFieldName = None):
    """
    Dot product
    """

    def ManipFunc(field, index):
      sum = 0.0
      for i, val in enumerate(field[index]):
        sum += val * array[index][i]
        
      return sum

    self.ManipulateField(fieldName, ManipFunc, newFieldName)

    return

  def CrossFieldWithField(self, fieldName, array, newFieldName = None, postMultiply = True):
    """
    Cross product
    """

    def ManipFunc(field, index):
      if postMultiply:
        return numpy.cross(field[index], array[index])
      else:
        return numpy.cross(array[index], field[index])

    self.ManipulateField(fieldName, ManipFunc, newFieldName)

    return

  def MatMulFieldWithField(self, fieldName, array, newFieldName = None, postMultiply = True):
    """
    Matrix multiplication
    """

    def ManipFunc(field, index):
      if postMultiply:
        return numpy.matrix(field[i]) * numpy.matix(array[i])
      else:
        return numpy.matix(array[i]) * numpy.matrix(field[i])

    self.ManipulateField(fieldName, ManipFunc, newFieldName)

    return

  # Default multiplication is dot product
  MulFieldByField = DotFieldWithField
  
  def GetDerivative(self, name):
    """
    Returns the derivative of field 'name', a
    vector field if 'name' is scalar, and a tensor field
    if 'name' is a vector. The field 'name' has to be point-wise data.
    The returned array gives a cell-wise derivative.
    """
    cd=vtk.vtkCellDerivatives()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      cd.SetInput(self.ugrid)
    else:
      cd.SetInputData(self.ugrid)
    pointdata=self.ugrid.GetPointData()
    nc=pointdata.GetArray(name).GetNumberOfComponents()
    if nc==1:
      cd.SetVectorModeToComputeGradient()
      cd.SetTensorModeToPassTensors()
      pointdata.SetActiveScalars(name)
      cd.Update()
      vtkdata=cd.GetUnstructuredGridOutput().GetCellData().GetArray('ScalarGradient')
      return arr([vtkdata.GetTuple3(i) for i in range(vtkdata.GetNumberOfTuples())])
    else:
      cd.SetTensorModeToComputeGradient()
      cd.SetVectorModeToPassVectors()
      pointdata.SetActiveVectors(name)
      cd.Update()
      vtkdata=cd.GetUnstructuredGridOutput().GetCellData().GetArray('VectorGradient')
      return arr([vtkdata.GetTuple9(i) for i in range(vtkdata.GetNumberOfTuples())])

  def GetVorticity(self, name):
    """
    Returns the vorticity of vectorfield 'name'.
    The field 'name' has to be point-wise data.
    The returned array gives a cell-wise derivative.
    """
    cd=vtk.vtkCellDerivatives()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      cd.SetInput(self.ugrid)
    else:
      cd.SetInputData(self.ugrid)
    pointdata=self.ugrid.GetPointData()
    cd.SetVectorModeToComputeVorticity()
    cd.SetTensorModeToPassTensors()
    pointdata.SetActiveVectors(name)
    cd.Update()
    vtkdata=cd.GetUnstructuredGridOutput().GetCellData().GetArray('VectorGradient')
    return arr([vtkdata.GetTuple3(i) for i in range(vtkdata.GetNumberOfTuples())])

  def CellDataToPointData(self):
    """
    Transforms all cell-wise fields in the vtu to point-wise fields.
    All existing fields will remain.
    """
    cdtpd=vtk.vtkCellDataToPointData()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      cdtpd.SetInput(self.ugrid)
    else:
      cdtpd.SetInputData(self.ugrid)
    cdtpd.PassCellDataOn()
    cdtpd.Update()
    self.ugrid=cdtpd.GetUnstructuredGridOutput()

class VTUProbe(object):
  """A class that combines a vtkProbeFilter with a list of invalid points (points that it failed to probe
  where we take the value of the nearest point)"""

  def __init__(self, ugrid, coordinates):
    # Initialise locator
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(ugrid)
    locator.SetTolerance(10.0)
    locator.Update()

    # Initialise probe
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    ilen, jlen = coordinates.shape
    for i in range(ilen):
      points.InsertNextPoint(coordinates[i][0], coordinates[i][1], coordinates[i][2])
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    self.probe = vtk.vtkProbeFilter()
    if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
      self.probe.SetInput(polydata)
      self.probe.SetSource(ugrid)
    else:
      self.probe.SetInputData(polydata)
      self.probe.SetSourceData(ugrid)
    self.probe.Update()

    # Generate a list invalidNodes, containing a map from invalid nodes in the
    # result to their closest nodes in the input
    valid_ids = self.probe.GetValidPoints()
    valid_loc = 0
    self.invalidNodes = []
    for i in range(ilen):
      if valid_ids.GetTuple1(valid_loc) == i:
        valid_loc += 1
      else:
        nearest = locator.FindClosestPoint([coordinates[i][0], coordinates[i][1], coordinates[i][2]])
        self.invalidNodes.append((i, nearest))
    self.ugrid = ugrid

  def GetField(self, name):
    # Get final updated values
    pointdata = self.probe.GetOutput().GetPointData()
    vtkdata=pointdata.GetArray(name)
    nc=vtkdata.GetNumberOfComponents()
    nt=vtkdata.GetNumberOfTuples()
    array = arr([vtkdata.GetValue(i) for i in range(nt * nc)])
    
    # Fix the point data at invalid nodes
    if len(self.invalidNodes) > 0:
      oldField = self.ugrid.GetPointData().GetArray(name)
      if oldField is None:
        oldField = self.ugrid.GetCellData().GetArray(name)
        if oldField is None:
          raise Exception("ERROR: couldn't find point or cell field data with name "+name+".")
      components = oldField.GetNumberOfComponents()
      for invalidNode, nearest in self.invalidNodes:
        for comp in range(nc):
          array[invalidNode * nc + comp] = oldField.GetValue(nearest * nc + comp)
          
    # this is a copy and paster from vtu.GetField above:
    if nc==9:
      return array.reshape(nt,3,3)
    elif nc==4:
      return array.reshape(nt,2,2)
    else:
      return array.reshape(nt,nc)
          
    return array
    
def VtuMatchLocations(vtu1, vtu2, tolerance = 1.0e-6):
  """
  Check that the locations in the supplied vtus match exactly, returning True if they
  match and False otherwise.
  The locations must be in the same order.
  """

  locations1 = vtu1.GetLocations().tolist()
  locations2 = vtu2.GetLocations()
  if not len(locations1) == len(locations2):
    return False
  for i in range(len(locations1)):
    if not len(locations1[i]) == len(locations2[i]):
      return False
    for j in range(len(locations1[i])):
      if abs(locations1[i][j] - locations2[i][j]) > tolerance:
        return False

  return True

def VtuMatchLocationsArbitrary(vtu1, vtu2, tolerance = 1.0e-6):
  """
  Check that the locations in the supplied vtus match, returning True if they
  match and False otherwise.
  The locations may be in a different order.
  """
   
  locations1 = vtu1.GetLocations()
  locations2 = vtu2.GetLocations()
  if not locations1.shape == locations2.shape:
    return False   
    
  for j in range(locations1.shape[1]):
    # compute the smallest possible precision given the range of this coordinate
    epsilon = numpy.finfo(numpy.float).eps * numpy.abs(locations1[:,j]).max()
    if tolerance<epsilon:
      # the specified tolerance is smaller than possible machine precision
      # (or something else went wrong)
      raise Exception("ERROR: specified tolerance is smaller than machine precision of given locations")
    # ensure epsilon doesn't get too small (might be for zero for instance)
    epsilon=max(epsilon,tolerance/100.0)

    # round to that many decimal places (-2 to be sure) so that
    # we don't get rounding issues with lexsort
    locations1[:,j]=numpy.around(locations1[:,j], int(-numpy.log10(epsilon))-2)
    locations2[:,j]=numpy.around(locations2[:,j], int(-numpy.log10(epsilon))-2)

  # lexical sort on x,y and z coordinates resp. of locations1 and locations2
  sort_index1=numpy.lexsort(locations1.T)
  sort_index2=numpy.lexsort(locations2.T)
  
  # should now be in same order, so we can check for its biggest difference
  return numpy.allclose(locations1[sort_index1],locations2[sort_index2], atol=tolerance)

def VtuDiff(vtu1, vtu2, fieldnamemap={}, original_fields=False):
  """
  Generate a vtu with fields generated by taking the difference between the field
  values in the two supplied vtus. Fields that are not common between the two vtus
  are neglected. If probe is True, the fields of vtu2 are projected onto the cell
  points of vtu1. Otherwise, the cell points of vtu1 and vtu2 must match.
  """

  # Generate empty output vtu
  resultVtu = vtu()

  # If the input vtu point locations match, do not use probe
  useProbe  = not VtuMatchLocations(vtu1, vtu2)
  if useProbe:
    probe = VTUProbe(vtu2.ugrid, vtu1.GetLocations())

  # Copy the grid from the first input vtu into the output vtu
  resultVtu.ugrid.DeepCopy(vtu1.ugrid)

  # Find common field names between the input vtus and generate corresponding
  # difference fields
  vtkdata=vtu1.ugrid.GetPointData()
  fieldnames1 = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
  vtkdata=vtu2.ugrid.GetPointData()
  fieldnames2 = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
  for fieldname1 in fieldnames1:
    fieldname2 = fieldnamemap.get(fieldname1, fieldname1)
    field1 = vtu1.GetField(fieldname1)
    if fieldname2 in fieldnames2:
      if useProbe:
        field2 = probe.GetField(fieldname2)
      else:
        field2 = vtu2.GetField(fieldname2)
      resultVtu.AddField(fieldname1, field1-field2)
      if original_fields: resultVtu.AddField(fieldname1+"::Original2", field2)
    else:
      resultVtu.RemoveField(fieldname1)
    if original_fields: resultVtu.AddField(fieldname1+"::Original1", field1)

  # Also look for cell-based fields. This only works if we don't have
  # to interpolate (both meshes are the same)
  vtkdata=vtu1.ugrid.GetCellData()
  fieldnames1 = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
  vtkdata=vtu2.ugrid.GetCellData()
  fieldnames2 = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
  if useProbe:
    # meshes are different - we can't interpolate cell-based fields so let's just remove them from the output
    for fieldname1 in fieldnames1:
      fieldname2 = fieldnamemap.get(fieldname1, fieldname1)
      if fieldname2=='vtkGhostLevels':
        # this field should just be passed on unchanged
        continue
      if original_fields: 
        field1 = vtu1.GetField(fieldname1)
        resultVtu.AddField(fieldname1+"::Original1", field1)
      resultVtu.RemoveField(fieldname1)
  else:
    # meshes are the same - we can simply subtract
    for fieldname1 in fieldnames1:
      fieldname2 = fieldnamemap.get(fieldname1, fieldname1)
      field1 = vtu1.GetField(fieldname1)
      if fieldname2=='vtkGhostLevels':
        # this field should just be passed on unchanged
        continue
      elif fieldname2 in fieldnames2:
        field2 = vtu2.GetField(fieldname2)
        resultVtu.AddField(fieldname1, field1-field2)
        if original_fields: resultVtu.AddField(fieldname1+"::Original2", field2)
      else:
        resultVtu.RemoveField(fieldname1)
      if original_fields: resultVtu.AddField(fieldname1+"::Original1", field1)

  return resultVtu

def extract_subvtus(filenames, meshfilename, region_ids, \
                    names=None, tindices=None, times=None, \
                    family=None, degree=None, writefile=False):
  """
  Extract vtus representing just the region_ids in region_ids.

  All other arguments are to deal with xdmfs and pvds with multiple times and/or functionspaces.
  """
  from buckettools import mesh, vtktools
  from lxml import etree
  import re 

  # sort out the files we're processing
  ltimes = []
  lfilenames = []
  for filename in filenames:
    basename, ext = os.path.splitext(filename)
    if ext == ".pvd":
      pvdtree = etree.parse(filename)
      pvdeles = [ele for ele in pvdtree.getroot().iterdescendants(tag="DataSet")]
      pvdtimes = numpy.asarray([float(ele.attrib["timestep"]) for ele in pvdeles])
      pvdfnames = [ele.attrib["file"] for ele in pvdeles]
      ltindices = list(range(len(pvdfnames)))
      if tindices is not None: ltindices = tindices.copy()
      if times is not None: ltindices = [numpy.abs(pvdtimes-t).argmin() for t in times]
      for ti in ltindices:
        ltimes.append(pvdtimes[ti])
        lfilenames.append(pvdfnames[ti])
    elif ext == ".xdmf":
      from buckettools import xdmftools
      # FIXME: this will write a whole load of vtus to disk
      xdmf = xdmftools.XDMF(filename)
      xdmffnames, xdmftimes = xdmf.vtus(names=names, tindices=tindices, times=times, \
                                        family=family, degree=degree, writefile=True)
      ltimes += xdmftimes
      lfilenames += xdmffnames
    elif ext in [".pvtu", ".vtu"]:
      try:
        vtutime = float(re.sub("[^0-9]", "", basename))
      except ValueError:
        vtutime = 0.0
      ltimes.append(vtutime)
      lfilenames.append(filename)
    else:
      raise Exception("ERROR: Unknown input file extension {} in extract_vtus!".format(ext,))

  # work out (default) output basename
  if isinstance(writefile, str):
    outbasename, outext = os.path.splitext(writefile)
    if outext not in ['', '.vtu', '.pvd']:
      raise Exception("ERROR: Unknown output file extension {} in extract_vtus!".format(ext,))

  # all of the files in lfilenames should now be [p]vtus

  # extract the submesh from the mesh described in meshfilename
  submesh, smcells, smfacets = mesh.extract_submesh(meshfilename, region_ids)
  cells = submesh.data().array("parent_cell_indices", 2)

  vtu = vtktools.vtu(lfilenames[0])
  # figure out which points are in the submesh cells in the vtu (not in the mesh!)
  points = []
  for ci in cells:
    cell = vtu.ugrid.GetCell(ci)
    np = cell.GetNumberOfPoints()
    for i in range(np):
      points.append(cell.GetPointId(i))
  points = list(set(points))

  # insert points into a new ugrid
  newugrid = vtktools.vtk.vtkUnstructuredGrid()
  newvtupoints = vtktools.vtk.vtkPoints()
  newvtupoints.SetDataTypeToDouble()
  newpoints = {}
  for i, pi in enumerate(points):
    newvtupoints.InsertNextPoint(vtu.ugrid.GetPoint(pi))
    newpoints[pi] = i
  newugrid.SetPoints(newvtupoints)

  # insert cells into a new grid
  for ci in cells:
    cell = vtu.ugrid.GetCell(ci)
    np = cell.GetNumberOfPoints()
    idList = vtktools.vtk.vtkIdList()
    for i in range(np): idList.InsertNextId(newpoints[cell.GetPointId(i)])
    newugrid.InsertNextCell(cell.GetCellType(), idList)

  # start looping over files and fields
  fi = 0
  newvtus = []
  for filename in lfilenames:
    vtu = vtktools.vtu(filename)

    newvtu = vtktools.vtu()
    newvtu.ugrid = newugrid

    fieldnames = vtu.GetFieldNames()
    for fname in fieldnames:
      field = vtu.GetField(fname)
      if field.shape[0] == vtu.ugrid.GetNumberOfPoints():
        newfield = [field[pi] for pi in points]
      elif field.shape[0] == vtu.ugrid.GetNumberOfCells():
        newfield = [field[i] for i in range(len(cells))]
      else:
        raise Exception("Unknown type of field!")
      newvtu.AddField(fname, numpy.asarray(newfield))

    if isinstance(writefile, str):
      loutbasename = outbasename
      if len(lfilenames) > 1: loutbasename += repr(fi).zfill(6)
    else:
      loutbasename = os.path.splitext(filename)[0]
      index = re.sub("[^0-9]", "", loutbasename)
      if len(index)>0: loutbasename = loutbasename[:-len(index)]
      loutbasename += "_subvtu"+index
    outfilename = loutbasename + ".vtu"
    if writefile is not False:
      newvtu.Write(outfilename)
      newvtus.append(outfilename)
    else:
      newvtu.filename = outfilename
      newvtus.append(newvtu)

    fi += 1

  pvdcollection = etree.Element('Collection')
  for ti, time in enumerate(ltimes):
    if isinstance(newvtus[ti], str):
      vtufilename = newvtus[ti]
    else:
      vtufilename = newvtus[ti].filename
    pvdcollection.append(etree.Element('DataSet', attrib={'timestep':repr(time), 'part':'0', 'file':vtufilename}))
  pvdroot = etree.Element('VTKFile', attrib={'type':'Collection', 'version':'0.1'})
  pvdroot.append(pvdcollection)
  pvdtree = etree.ElementTree(pvdroot)

  if writefile is not False:
    if isinstance(writefile, str): 
      loutbasename = outbasename
    else:
      loutbasename = os.path.splitext(lfilenames[0])[0]
      index = re.sub("[^0-9]", "", loutbasename)
      loutbasename = loutbasename[:-len(index)]+"_subvtus"
    pvdtree.write(loutbasename+".pvd", encoding='UTF-8', pretty_print=True, xml_declaration=True)
    pvd = basename+".pvd"
  else:
    pvd = pvdtree

  return pvd, newvtus, ltimes

def timeaverage_vtu(filenames, \
                    names=None, tindices=None, times=None, \
                    family=None, degree=None, writefile=False):
  """
  Generate a time-averaged vtu of the files in filenames.

  All other arguments are to deal with xdmfs and pvds with multiple times and/or functionspaces.
  """
  from buckettools import mesh, vtktools
  from lxml import etree
  import re 

  # sort out the files we're processing
  ltimes = []
  lfilenames = []
  for filename in filenames:
    basename, ext = os.path.splitext(filename)
    if ext == ".pvd":
      pvdtree = etree.parse(filename)
      pvdeles = [ele for ele in pvdtree.getroot().iterdescendants(tag="DataSet")]
      pvdtimes = numpy.asarray([float(ele.attrib["timestep"]) for ele in pvdeles])
      pvdfnames = [ele.attrib["file"] for ele in pvdeles]
      ltindices = list(range(len(pvdfnames)))
      if tindices is not None: ltindices = tindices.copy()
      if times is not None: ltindices = [numpy.abs(pvdtimes-t).argmin() for t in times]
      for ti in ltindices:
        ltimes.append(pvdtimes[ti])
        lfilenames.append(pvdfnames[ti])
    elif ext == ".xdmf":
      from buckettools import xdmftools
      # FIXME: this will write a whole load of vtus to disk
      xdmf = xdmftools.XDMF(filename)
      xdmffnames, xdmftimes = xdmf.vtus(names=names, tindices=tindices, times=times, \
                                        family=family, degree=degree, writefile=True)
      ltimes += xdmftimes
      lfilenames += xdmffnames
    elif ext in [".pvtu", ".vtu"]:
      vtutime = float(re.sub("[^0-9]", "", basename))
      ltimes.append(vtutime)
      lfilenames.append(filename)
    else:
      raise Exception("ERROR: Unknown input file extension {} in timeaverage_vtu!".format(ext,))

  # all of the files in lfilenames should now be [p]vtus

  vtu = vtktools.vtu(lfilenames[0])
  # figure out fields we'll be averaging
  fieldnames = vtu.GetFieldNames()
  fields = {}
  oldfields = {}
  for fname in fieldnames:
    field = vtu.GetField(fname)
    fields[fname] = numpy.zeros(field.shape)
    oldfields[fname] = field

  oldtime = ltimes[0]

  # start looping over files and fields
  for fi, filename in enumerate(lfilenames[1:]):
    try:
      vtu = vtktools.vtu(filename)
    except:
      continue

    time = ltimes[fi+1]
    dt = time-oldtime

    for fname, field in fields.items():
      field = vtu.GetField(fname)
      fields[fname] += 0.5*dt*(field + oldfields[fname])
      oldfields[fname] = field

    oldtime = time

  vtu = vtktools.vtu(lfilenames[0])

  # figure out fields we'll be averaging
  for fname, field in fields.items(): vtu.AddField(fname, field)
  
  # output if requested
  if writefile is not False:
    if isinstance(writefile, str):
      outbasename, outext = os.path.splitext(writefile)
      if outext not in ['', '.vtu']:
        raise Exception("ERROR: Unknown output file extension {} in timeaverage_vtu!".format(ext,))
    else:
      outbasename = os.path.splitext(lfilenames[0])[0]
      index = re.sub("[^0-9]", "", outbasename)
      outbasename = outbasename[:-len(index)]+"_timeaveraged"
    vtu.Write(outbasename+".vtu")
    return outbasename+".vtu"
  else:
    return vtu

