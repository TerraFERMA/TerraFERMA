#!/usr/bin/python

import dolfin
from optparse import OptionParser
import sys
import numpy
import pylab

optparser=OptionParser(usage='usage: %prog <filename>',
                       add_help_option=True,
                       description="""hmmm""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

filename = argv[0]

dolfin.set_log_level(1)

# process the dolfin (xml) mesh file
xml_mesh = dolfin.Mesh(filename+".xml")
xml_mesh_test = dolfin.Mesh(filename+".xml")
xml_mesh.init()
assert((xml_mesh_test.coordinates() == xml_mesh.coordinates()).all())
assert((xml_mesh_test.cells()==xml_mesh.cells()).all())

cell_subdomain = dolfin.MeshFunction("uint", xml_mesh, xml_mesh.topology().dim())
cell_subdomain.set_all(0)
edge_subdomain = dolfin.MeshFunction("uint", xml_mesh, xml_mesh.topology().dim()-1)
edge_subdomain.set_all(0)

topo = xml_mesh.topology()
edge2vertex = topo(1,0)
cell2vertex = topo(2,0)

# process the gmsh mesh file
msh_mesh = file(filename+".msh", 'r')
# Header section
assert(msh_mesh.readline().strip()=="$MeshFormat")
assert(msh_mesh.readline().strip()in["2 0 8", "2.1 0 8", "2.2 0 8"])
assert(msh_mesh.readline().strip()=="$EndMeshFormat")

# Nodes section
while msh_mesh.readline().strip() !="$Nodes":
   pass

msh_num_vertices=int(msh_mesh.readline())
assert(msh_num_vertices==xml_mesh.num_vertices())

# Skip the nodes
while msh_mesh.readline().strip() !="$EndNodes":
   pass

# Elements section
assert(msh_mesh.readline().strip()=="$Elements")
msh_num_elements = int(msh_mesh.readline())
# note that this won't necessarily be the same as xml_mesh.num_cells() because
# it hopefully includes the surface information we're looking to retrieve

# Now loop over the elements placing them in the appropriate buckets.
msh_edges=[]
msh_triangles=[]
msh_tets=[]
#msh_quads=[]
#msh_hexes=[]

for i in range(msh_num_elements):
    element=msh_mesh.readline().split()
    if (element[1]=="1"):
        msh_edges.append([int(e)-1 for e in element[-2:]]+[int(element[3])])
    elif (element[1]=="2"):
        msh_triangles.append([int(e)-1 for e in element[-3:]]+[int(element[3])])
#    elif (element[1]=="3"):
#        msh_quads.append([int(e)-1 for e in element[-4:]]+[int(element[3])])
    elif (element[1]=="4"):
        msh_tets.append([int(e)-1 for e in element[-4:]]+[int(element[3])])
#    elif (element[1]=="5"):
#        msh_hexes.append([int(e)-1 for e in element[-8:]]+[int(element[3])])
    elif(element[1]=="15"):
        # Ignore point elements
        pass
    else:
        sys.stderr.write("Unknown element type "+`element[1]`+'\n')
        sys.exit(1)

#if len(msh_tets) > 0:
#  if len(msh_hexes) > 0:
#    sys.stderr.write("Warning: Mixed tet/hex mesh encountered - discarding hexes")
#  if len(msh_quads) > 0:
#    sys.stderr.write("Warning: Mixed tet/quad mesh encountered - discarding quads")
#elif len(msh_triangles) > 0:
#  if len(msh_hexes) > 0:
#    sys.stderr.write("Warning: Mixed triangle/hex mesh encountered - discarding hexes")
#  if len(msh_quads) > 0:
#    sys.stderr.write("Warning: Mixed triangle/quad mesh encountered - discarding quads")

if(edge_subdomain.dim()==1):
  msh_facets = numpy.array(msh_edges)
elif(edge_subdomain.dim()==2):
  msh_facets = numpy.array(msh_triangles)
else:
  sys.stderr.write("My universe is 3D, sorry Albert.")
  sys.exit(1)

if(cell_subdomain.dim()==2):
  msh_cells = numpy.array(msh_triangles)
elif(cell_subdomain.dim()==3):
  msh_cells = numpy.array(msh_tets)
else:
  sys.stderr.write("Can only cope with 2 or 3D meshes, sorry.")
  sys.exit(1)
assert(msh_cells.shape[0]==xml_mesh.num_cells())

# do a lot of searching for edges
for edge in range(xml_mesh.num_edges()):
  indices = numpy.array([edge2vertex(edge).size*i for i in range(msh_facets.shape[0])])
  for i in range(edge2vertex(edge).size):
    indices_i = pylab.find(edge2vertex(edge)[i]==msh_facets[indices/edge2vertex(edge).size,:-1])
    indices = indices[indices_i/edge2vertex(edge).size]
    if(indices.size==0): break # no point in still searching if we failed to find any indices
  if ((i==edge2vertex(edge).size-1) and (indices.size>0)):
    # we made it to the end of the loop and we found an index
    if (indices.size==1):
      # only one index mind
      edge_subdomain[edge] = msh_facets[indices[0]/edge2vertex(edge).size, -1]
    else:
      # or perhaps more... oh dear
      sys.stderr.write("Found more than one matching facet.")
      sys.exit(1)

edge_file = dolfin.File(filename+"_edge_subdomain.xml")
edge_file << edge_subdomain

# do a lot of searching for cells
for cell in range(xml_mesh.num_cells()):
  indices = numpy.array([cell2vertex(cell).size*i for i in range(msh_cells.shape[0])])
  for i in range(cell2vertex(cell).size):
    indices_i = pylab.find(cell2vertex(cell)[i]==msh_cells[indices/cell2vertex(cell).size,:-1])
    indices = indices[indices_i/cell2vertex(cell).size]
    if(indices.size==0): break # no point in still searching if we failed to find any indices
  if ((i==cell2vertex(cell).size-1) and (indices.size>0)):
    # we made it to the end of the loop and we found an index
    if (indices.size==1):
      # only one index mind
      cell_subdomain[cell] = msh_cells[indices[0]/cell2vertex(cell).size, -1]
    else:
      # or perhaps more... oh dear
      sys.stderr.write("Found more than one matching cell.")
      sys.exit(1)

cell_file = dolfin.File(filename+"_cell_subdomain.xml")
cell_file << cell_subdomain

