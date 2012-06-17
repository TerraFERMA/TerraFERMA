
import numpy
from scipy import interpolate as interp
from scipy import optimize as opt
from math import sqrt
import pylab
import sys

class GeoFile:
  def __init__(self):
    self.pindex = 1
    self.cindex = 1
    self.sindex = 1
    self.lines = []
  
  def addpoint(self, point, comment=None):
    if point.eid:
      return
    point.eid = self.pindex
    self.pindex += 1
    line = "Point("+`point.eid`+") = {"+`point.x`+", "+`point.y`+", "+`point.z`+", "+`point.res`+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += "\n"
    self.lines.append(line)

  def addcurve(self, curve, comment=None):
    if curve.eid:
      return
    curve.eid = self.cindex
    self.cindex += 1
    for p in curve.points:
      self.addpoint(p)
    line = curve.type+"("+`curve.eid`+") = {"
    for p in range(len(curve.points)-1):
      line += `curve.points[p].eid`+", "
    line += `curve.points[-1].eid`+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += "\n"
    self.lines.append(line)

  def addsurface(self, surface, comment=None):
    if surface.eid:
      return
    surface.eid = self.sindex
    self.sindex += 1
    for curve in surface.curves:
      self.addcurve(curve)
    line = "Line Loop("+`self.cindex`+") = {"
    for c in range(len(surface.curves)-1):
      line += `surface.directions[c]*surface.curves[c].eid`+", "
    line += `surface.directions[-1]*surface.curves[-1].eid`+"};\n"
    self.lines.append(line)
    line = "Plane Surface("+`surface.eid`+") = {"+`self.cindex`+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += "\n"
    self.lines.append(line)
    self.cindex += 1

  def addphysicalpoint(self, pid, points):
    line = "Physical Point("+`pid`+") = {"
    for p in range(len(points)-1):
      line += `points[p].eid`+", "
    line += `points[-1].eid`+"};\n"
    self.lines.append(line)

  def addphysicalline(self, pid, curves):
    line = "Physical Line("+`pid`+") = {"
    for c in range(len(curves)-1):
      line += `curves[c].eid`+", "
    line += `curves[-1].eid`+"};\n"
    self.lines.append(line)

  def addphysicalsurface(self, pid, surfaces):
    line = "Physical Surface("+`pid`+") = {"
    for s in range(len(surfaces)-1):
      line += `surfaces[s].eid`+", "
    line += `surfaces[-1].eid`+"};\n"
    self.lines.append(line)

  def addembed(self, surface, items):
    line = items[0].type+" {"
    for i in range(len(items)-1):
      assert(items[i].type==items[0].type)
      line += `items[i].eid`+", "
    line += `items[-1].eid`+"} In Surface {"+`surface.eid`+"};\n"
    self.lines.append(line) 

  def linebreak(self):
    self.lines.append("\n")

  def comment(self, string):
    return "// "+string+"\n"

  def produced_comment(self):
    self.lines.append(self.comment("Produced by: "+" ".join(sys.argv)))

  def write(self, filename):
    f = file(filename, 'w')
    f.writelines(self.lines)
    f.close()

class ElementaryEntity:
  def __init__(self, name=None):
    self.name = name
    self.eid  = None

class PhysicalEntity:
  def __init__(self, pid=None):
    self.pid = pid

class Point(ElementaryEntity, PhysicalEntity):
  def __init__(self, coord, name=None, res=1.0):
    self.type = "Point"
    self.x = coord[0]
    if len(coord) > 1: 
      self.y = coord[1]
    else:
      self.y = 0.0
    if len(coord) > 2: 
      self.z = coord[2]
    else:
      self.z = 0.0
    self.res = res
    ElementaryEntity.__init__(self, name=name)
    PhysicalEntity.__init__(self)

class Curve(ElementaryEntity, PhysicalEntity):
  def __init__(self, points, name=None, pid=None):
    self.points = points
    self.name = name
    self.x = None
    self.y = None
    self.u = None
    ElementaryEntity.__init__(self, name=name)
    PhysicalEntity.__init__(self, pid=pid)

  def update(self):
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    isort = numpy.argsort(self.x)
    points = numpy.asarray(self.points)[isort]
    self.points = points.tolist()
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    self.y = numpy.array([self.points[i].y for i in range(len(self.points))])

class Line(Curve):
  def __init__(self, points, name=None, pid=None):
    assert(len(points)==2)
    self.type = "Line"
    Curve.__init__(self, points, name=name, pid=pid)
    self.update()

  def update(self):
    Curve.update(self)
    self.u = numpy.arange(0.0, self.x.size)

  def __call__(self, u):
    return [self.points[0].x + u*(self.points[1].x - self.points[0].x), \
            self.points[0].y + u*(self.points[1].y - self.points[0].y)]

class Spline(Curve):
  def __init__(self, points, name=None, pid=None):
    self.type = "Spline"
    Curve.__init__(self, points, name=name, pid=pid)
    self.cmr = None
    self.update()

  def __call__(self, u, der=0):
    return [self.cmr[0].derivative(u,der=der), self.cmr[1].derivative(u,der=der)]

  def update(self):
    Curve.update(self)
    self.cmr, self.u = self.catmullrom()

  def find_derivatives(self):
    dk = [numpy.zeros_like(self.x) for i in range(2)]
    dk[0][0]    =  self.x[1]  - self.x[0]
    dk[0][1:-1] = (self.x[2:] - self.x[:-2])/2.0
    dk[0][-1]   =  self.x[-1] - self.x[-2]
    dk[1][0]    =  self.y[1]  - self.y[0]
    dk[1][1:-1] = (self.y[2:] - self.y[:-2])/2.0
    dk[1][-1]   =  self.y[-1] - self.y[-2]
    return dk

  def catmullrom(self):
    derivs = self.find_derivatives()
    u = numpy.arange(0.0, self.x.size)
    polys = []
    polys.append(interp.PiecewisePolynomial(u, zip(self.x, derivs[0]), orders=3, direction=None))
    polys.append(interp.PiecewisePolynomial(u, zip(self.y, derivs[1]), orders=3, direction=None))
    return polys,u

  def intersecty(self, yint):
    loc = abs(self.y-yint).argmin()
    if loc == 0:
      loc0 = loc
      loc1 = loc+1
    elif loc == self.y.size-1:
      loc0 = loc-1
      loc1 = loc
    else:
      if (self.y[loc]-yint)*(self.y[loc+1]-yint) < 0.0:
        loc0 = loc
        loc1 = loc+1
      else:
        loc0 = loc-1
        loc1 = loc
    uint = opt.brentq(lambda p: self.cmr[1](p)-yint, self.u[loc0], self.u[loc1])
    return self(uint)
      
  def intersectx(self, xint):
    loc = abs(self.x-xint).argmin()
    if loc == 0:
      loc0 = loc
      loc1 = loc+1
    elif loc == self.x.size-1:
      loc0 = loc-1
      loc1 = loc
    else:
      if (self.x[loc]-xint)*(self.x[loc+1]-xint) < 0.0:
        loc0 = loc
        loc1 = loc+1
      else:
        loc0 = loc-1
        loc1 = loc
    uint = opt.brentq(lambda p: self.cmr[0](p)-xint, self.u[loc0], self.u[loc1])
    return self(uint)
    
  def translatenormal(self, dist):
    for i in range(len(self.u)):
      der = self(self.u[i], der=1)
      vec = numpy.array([-der[1], der[0]])
      vec = vec/sqrt(sum(vec**2))
      self.points[i].x += dist*vec[0]
      self.points[i].y += dist*vec[1]
    self.update()

  def croppoint(self, pind, coord):
    for i in range(len(pind)-1):
      p = self.points.pop(pind[i])
    self.points[pind[-1]].x = coord[0]
    self.points[pind[-1]].y = coord[1]
    self.update()

  def crop(self, left=None, bottom=None, right=None, top=None):
    if left:
      out = pylab.find(self.x < left)
      if (len(out)>0):
        coord = self.intersectx(left)
        self.croppoint(out, coord)
    if right:
      out = pylab.find(self.x > right)
      if (len(out)>0):
        coord = self.intersectx(right)
        self.croppoint(out, coord)
    if bottom:
      out = pylab.find(self.y < bottom)
      if (len(out)>0):
        coord = self.intersecty(bottom)
        self.croppoint(out, coord)
    if top:
      out = pylab.find(self.y > top)
      if (len(out)>0):
        coord = self.intersecty(top)
        self.croppoint(out, coord)

  def translatenormalandcrop(self, dist):
    left = self.x.min()
    right = self.x.max()
    bottom = self.y.min()
    top = self.y.max()
    self.translatenormal(dist)
    self.crop(left=left, bottom=bottom, right=right, top=top)

  def findpointx(self, xint):
    ind = pylab.find(self.x==xint)
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpointy(self, yint):
    ind = pylab.find(self.y==yint)
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpoint(self, name):
    for p in self.points:
      if p.name == name: return p
    return None

  def appendpoint(self, p):
    self.points.append(p)
    self.update()

  def split(self, ps):
    points = [point for point in self.points]
    pind = range(len(points))
    for p in ps:
      if p not in points or p==points[0] or p==points[-1]: continue
      for i in range(1,len(points)-1):
        if points[i] == p: break
      newp0 = Point(self(self.u[pind[i]]-0.1), res=p.res)
      newp1 = Point(self(self.u[pind[i]]+0.1), res=p.res)
      points.insert(i, newp0)
      pind.insert(i, None)
      points.insert(i+2, newp1)
      pind.insert(i+2, None)
    splines = []
    j = 0
    for p in range(len(ps)+1):
      newpoints = []
      for i in range(j,len(points)):
        newpoints.append(points[i])
        if p<len(ps):
          if points[i] == ps[p]: break
      j = i
      if self.name:
        name = self.name+"_split"+`p`
      else: 
        name = None
      splines.append(Spline(newpoints, name=name, pid=self.pid))
    return splines

class Surface(ElementaryEntity, PhysicalEntity):
  def __init__(self, curves, name=None, pid = None):
    self.type = "Surface"
    self.curves = [curves[0]]
    self.directions = [1]
    ind = range(1, len(curves))
    for j in range(1, len(curves)):
      pcurve = self.curves[j-1]
      if self.directions[j-1]==-1:
        pind = 0
      else:
        pind = -1
      for i in range(len(ind)):
        if pcurve.points[pind] == curves[ind[i]].points[0]:
          self.directions.append(1)
          break
        elif pcurve.points[pind] == curves[ind[i]].points[-1]:
          self.directions.append(-1)
          break
      self.curves.append(curves[ind[i]])
      i = ind.pop(i)
    ElementaryEntity.__init__(self, name=name)
    PhysicalEntity.__init__(self, pid=pid)

class Geometry:
  def __init__(self):
    self.namedpoints = {}
    self.namedcurves = {}
    self.namedsurfaces = {}
    self.physicalpoints = {}
    self.physicalcurves = {}
    self.physicalsurfaces = {}
    self.points  = []
    self.curves  = []
    self.surfaces = []
    self.pointcomments  = []
    self.curvecomments  = []
    self.surfacecomments = []
    self.pointembeds = {}
    self.lineembeds = {}
    self.geofile = None

  def addpoint(self, point, name=None, comment=None):
    self.points.append(point)
    self.pointcomments.append(comment)
    if not comment:
      if name:
        self.pointcomments[-1] = name
      elif point.name:
        self.pointcomments[-1] = point.name
    if name:
      self.namedpoints[name] = point
    elif point.name:
      self.namedpoints[point.name] = point
    if point.pid:
      if point.pid in self.physicalpoints:
        self.physicalpoints[point.pid].append(point)
      else:
        self.physicalpoints[point.pid] = [point]

  def addcurve(self, curve, name=None, comment=None):
    self.curves.append(curve)
    self.curvecomments.append(comment)
    if not comment:
      if name:
        self.curvecomments[-1] = name
      elif curve.name:
        self.curvecomments[-1] = curve.name
    if name:
      self.namedcurves[name] = curve
    elif curve.name:
      self.namedcurves[curve.name] = curve
    if curve.pid:
      if curve.pid in self.physicalcurves:
        self.physicalcurves[curve.pid].append(curve)
      else:
        self.physicalcurves[curve.pid] = [curve]

  def addsurface(self, surface, name=None, comment=None):
    self.surfaces.append(surface)
    self.surfacecomments.append(comment)
    if not comment:
      if name:
        self.surfacecomments[-1] = name
      elif surface.name:
        self.surfacecomments[-1] = surface.name
    if name:
      self.namedsurfaces[name] = surface
    elif surface.name:
      self.namedsurfaces[surface.name] = surface
    if surface.pid:
      if surface.pid in self.physicalsurfaces:
        self.physicalsurfaces[surface.pid].append(surface)
      else:
        self.physicalsurfaces[surface.pid] = [surface]

  def addembed(self, surface, item):
    assert(surface.type=="Surface")
    assert(item.type=="Point" or item.type=="Line")
    if item.type=="Point":
      if surface in self.pointembeds:
        self.pointembeds[surface].append(item)
      else:
        self.pointembeds[surface] = [item]
    else:
      if surface in self.lineembeds:
        self.lineembeds[surface].append(item)
      else:
        self.lineembeds[surface] = [item]

  def writegeofile(self, filename):
    for surface in self.surfaces:
      surface.eid = None
      for curve in surface.curves:
        curve.eid = None
        for point in curve.points:
          point.eid = None
    for curve in self.curves:
      curve.eid = None
      for point in curve.points:
        point.eid = None
    for point in self.points:
      point.eid = None

    self.geofile = GeoFile()
    for s in range(len(self.surfaces)):
      surface = self.surfaces[s]
      self.geofile.addsurface(surface, self.surfacecomments[s])
      self.geofile.linebreak()
    self.geofile.linebreak()
    for c in range(len(self.curves)):
      curve = self.curves[c]
      self.geofile.addcurve(curve, self.curvecomments[c])
      self.geofile.linebreak()
    self.geofile.linebreak()
    for p in range(len(self.points)):
      point = self.points[p]
      self.geofile.addpoint(point, self.pointcomments[p])
    self.geofile.linebreak()

    for surface,items in self.pointembeds.iteritems():
      self.geofile.addembed(surface,items)
    for surface,items in self.lineembeds.iteritems():
      self.geofile.addembed(surface,items)
    self.geofile.linebreak()
    
    for pid,surfaces in self.physicalsurfaces.iteritems():
      self.geofile.addphysicalsurface(pid,surfaces)
    self.geofile.linebreak()
    for pid,curves in self.physicalcurves.iteritems():
      self.geofile.addphysicalline(pid,curves)
    self.geofile.linebreak()
    for pid,points in self.physicalpoints.iteritems():
      self.geofile.addphysicalpoint(pid,points)
    self.geofile.linebreak()

    self.geofile.produced_comment()

    self.geofile.write(filename)

  def pylabplot(self, lineres=100):
    for curve in self.curves:
      unew = numpy.arange(curve.u[0], curve.u[-1]+((curve.u[-1]-curve.u[0])/(2.*lineres)), 1./lineres)
      pylab.plot(curve(unew)[0], curve(unew)[1], 'k')
      pylab.plot(curve.x, curve.y, 'ok')
    for point in self.points:
      pylab.plot(point.x, point.y, 'ok')
  
