import dolfin
import numpy
import subprocess

def makemesh(vertex_coordinates, cells):
    "constructs dolfin mesh given numpy arrays of vertex coordinates and cell maps, *does not* output mesh to file"

    # create mesh and editor
    mesh = dolfin.Mesh()
    M = dolfin.MeshEditor()

    # get information about mesh
    num_vertices, gdim = vertex_coordinates.shape
    num_cells, cdim = cells.shape
    tdim = cdim -  1 # topological dimension
    if tdim == 2:
        shape="triangle"
    elif tdim == 3:
        shape = "tetrahedron"
    else:
        print "Error: unknown cell shape"
        return

    # make sure cell vertex indices are 0 offset (C style) rather than 1 offset (Fortran style)
    cells -= cells.min()


    # initialize mesh
    M.open(mesh,shape,tdim,gdim)
    M.init_vertices(num_vertices)
    M.init_cells(num_cells)

    # add vertices
    for i in range(num_vertices):
        p = dolfin.Point(gdim,vertex_coordinates[i,:])
        M.add_vertex(i,p)

    # add cells
    for i in range(num_cells):
        #c = STLVectorUInt(cells[i,:])
        M.add_cell(i, cells[i,0], cells[i,1], cells[i,2])

    # finish up and order
    M.close()

    return mesh

def writemesh(filename, vertex_coordinates, cells):
    "constructs dolfin mesh given numpy arrays of vertex coordinates and cell maps, outputs mesh to file"

    # create mesh and editor
    mesh = makemesh(vertex_coordinates, cells)

    # write mesh (compressed or uncompressed)
    if filename.find('.xml.gz') > 0:
        name = filename.replace('.gz','')
        f=dolfin.File(name)
        f << mesh
        subprocess.call(['gzip','-vf',name])
    else:
        if filename.find('.xml'):
            f=dolfin.File(filename)
            f << mesh
        else:
            print "Error: %s does not appear to be an xml file" % filename

    return mesh

def splitmesh(mesh, split_facet_ids, centroid_function, preserve_facet_ids=[], new_facet_ids=None):
  # get the indices of facets (those to be split, those to be preserved and all indexed facets)
  facetsubdomain = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim()-1, mesh.domains())
  split_facet_indices = numpy.where([fid in split_facet_ids for fid in facetsubdomain.array()])[0]
  preserve_facet_indices = numpy.where([fid in preserve_facet_ids for fid in facetsubdomain.array()])[0]
  all_facet_indices = numpy.where(facetsubdomain.array()!=0)[0]

  # duplicate the cell data and the coordinate data
  split_cells = numpy.copy(mesh.cells())
  split_coords = numpy.copy(mesh.coordinates())

  # set up the topological mappings we'll need below
  topo = mesh.topology()
  facets2vertices = topo(topo.dim()-1,0)
  vertices2facets = topo(0,topo.dim()-1)
  vertices2cells  = topo(0,topo.dim())

  # set up a set of the vertex indices that need to be split
  new_coords_ind = set()
  # loop over the facet indices that we've requested to be split
  for i in range(len(split_facet_indices)):
    # loop over the vertices that make up this facet
    for ind in facets2vertices(split_facet_indices[i]):
      # check that the current vertex isn't attached to any facet that we've been told
      # to preserve
      if not numpy.any([f in preserve_facet_indices for f in vertices2facets(ind)]):
        # if not, record the index of the vertex
        new_coords_ind.add(ind)

  # set up a map from the index of the old vertex to the index of the new, duplicated vertex
  new_coords_map = []
  for ind in new_coords_ind:
    print "Splitting vertex at ", split_coords[ind]
    # duplicate the coordinates of the vertex
    split_coords = numpy.concatenate((split_coords, numpy.array([split_coords[ind]])), axis=0)
    # map from the old vertex to the new one
    new_coords_map.append([ind, len(split_coords)-1])

  # turn the list into a map
  new_coords_map = dict(new_coords_map)
  # and invert the map so that you can jump from the new vertex index to the old
  old_coords_map = dict((v,k) for k, v in new_coords_map.iteritems())

  # keep track of which cells we've visited
  visited = numpy.zeros(mesh.num_cells(), numpy.int)

  # loop over the vertices that are being split again
  for ind in new_coords_ind:
    # loop over the cells attached to these vertices
    for cell_ind in vertices2cells(ind):
      # if this cell hasn't been visited before
      if(not visited[cell_ind]):
        # mark it as having been visited
        visited[cell_ind] = 1
        # work out the cell centroid
        cell_centroid = mesh.coordinates()[mesh.cells()[cell_ind]].sum(axis=0)/float(mesh.cells()[cell_ind].size)
        # if the centroid is above the x==y line then switch this cell to using the new split vertices
        if(centroid_function(cell_centroid)):
          ## need to import pylab if these lines are uncommented
          # pylab.figure(0)
          # pylab.plot(numpy.append(mesh.coordinates()[mesh.cells()[cell_ind]][:,0], mesh.coordinates()[mesh.cells()[cell_ind]][0,0]),
          #            numpy.append(mesh.coordinates()[mesh.cells()[cell_ind]][:,1], mesh.coordinates()[mesh.cells()[cell_ind]][0,1]))
          # loop over the vertices attached to this cell (still the old ones)
          for i in range(len(mesh.cells()[cell_ind])):
            if mesh.cells()[cell_ind][i] in new_coords_ind: 
              # replace those vertices that are being split in the mesh
              split_cells[cell_ind][i] = new_coords_map[mesh.cells()[cell_ind][i]]
              assert(split_cells[cell_ind][i] != mesh.cells()[cell_ind][i])
          # pylab.figure(1)
          # pylab.plot(numpy.append(split_coords[split_cells[cell_ind]][:,0], split_coords[split_cells[cell_ind]][0,0]),
          #            numpy.append(split_coords[split_cells[cell_ind]][:,1], split_coords[split_cells[cell_ind]][0,1]))

  # pylab.show()

  # write the new split mesh to file
  split_mesh = makemesh(split_coords, split_cells)
  split_mesh.init()

  # initialize the topology of the split mesh
  split_topo = split_mesh.topology()
  split_vertices2facets = split_topo(0,split_topo.dim()-1)
  split_facets2vertices = split_topo(split_topo.dim()-1,0)

  # initialize the new mesh functions
  split_mesh.domains().init(split_topo.dim())
  cellmarkers = mesh.domains().markers(split_topo.dim())
  for k,v in cellmarkers.iteritems():
    split_mesh.domains().set_marker((k, v), split_topo.dim())

  # FIXME: this could be done better by just creating a MeshValueCollection as
  # we loop over cells in the previous iteration (wasn't available when this
  # function was first written!

  # loop over all facets with non-zero ids in the original mesh
  for facet in all_facet_indices:
    vertices = facets2vertices(facet)
    # loop over the vertices in the current facet
    for vertex in vertices:
      # loop over the facets connected to this vertex in the new, split mesh 
      for split_facet in split_vertices2facets(vertex):
        # loop over the vertices in the facet in the new, split mesh
        for split_neigh in split_facets2vertices(split_facet):
          # if the neighbouring vertex index is not the current base vertex (i.e. it is a neighbour)
          # and the neighbour is a vertex in the original facet
          if((split_neigh != vertex) and (split_neigh in vertices)):
            # then label this facet with the original facet id
            split_mesh.domains().set_marker((int(split_facet), facetsubdomain[int(facet)]), split_topo.dim()-1)
      # if the vertex is one of the ones that was split
      if vertex in new_coords_ind:
        # find out the new vertex number
        split_vertex = new_coords_map[vertex]
        # then sloop over the facets connected to the new, split vertex
        for split_facet in split_vertices2facets(split_vertex):
          # loop over the neighbouring vertices
          for split_neigh in split_facets2vertices(split_facet):
            # if this really is a neighbour (i.e. it isn't the original split vertex)
            if(split_neigh!=split_vertex):
              # if the neighbour was split itself
              if(split_neigh in old_coords_map):
                # find out the original index for the neighbour
                neigh = old_coords_map[split_neigh]
                # if it's one of the vertices in the current facet
                if (neigh in vertices):
                  # assign this facet a new facet id
                  assert(facetsubdomain[int(facet)] in split_facet_ids)
                  if new_facet_ids:
                    # user has supplied new ids
                    assert(len(numpy.where([fid==facetsubdomain[int(facet)] for fid in split_facet_ids])[0])==1)
                    newid = new_facet_ids[numpy.where([fid==facetsubdomain[int(facet)] for fid in split_facet_ids])[0][0]]
                  else:
                    # user hasn't supplied new ids
                    newid = int(`facetsubdomain[int(facet)]`+'00')
                  split_mesh.domains().set_marker((int(split_facet), newid), split_topo.dim()-1)
              # the neighbour wasn't itself split
              else:
                # so the neighbour has the same index
                neigh = split_neigh
                # but if it's in the vertices this means we're at the end of the split
                if (neigh in vertices):
                  print 'found final split facet, vertex, split_vertex, neigh ', vertex, split_vertex, neigh, \
                                                                                mesh.coordinates()[vertex], \
                                                                                mesh.coordinates()[neigh]
                  # if this was an id that we were splitting, give it a new facet id
                  if facetsubdomain[int(facet)] in split_facet_ids:
                    if new_facet_ids:
                      # user has supplied new ids
                      assert(len(numpy.where([fid==facetsubdomain[int(facet)] for fid in split_facet_ids])[0])==1)
                      newid = new_facet_ids[numpy.where([fid==facetsubdomain[int(facet)] for fid in split_facet_ids])[0][0]]
                    else:
                      # user hasn't supplied new ids
                      newid = int(`facetsubdomain[int(facet)]`+'00')
                    print '  providing new id', newid
                    split_mesh.domains().set_marker((int(split_facet), newid), split_topo.dim()-1)
                  # otherwise give it the same old facet id
                  else:
                    print '  using old id'
                    split_mesh.domains().set_marker((int(split_facet), facetsubdomain[int(facet)]), split_topo.dim()-1)

  return split_mesh 
