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


import numpy
from scipy import interpolate as interp
from scipy import integrate as integ
from scipy import optimize as opt
from math import sqrt
import sys
import os

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
    line = "Point("+repr(point.eid)+") = {"+repr(point.x)+", "+repr(point.y)+", "+repr(point.z)+", "+repr(point.res)+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += os.linesep
    self.lines.append(line)

  def addcurve(self, curve, comment=None):
    if curve.eid:
      return
    curve.eid = self.cindex
    self.cindex += 1
    for p in curve.points:
      self.addpoint(p)
    line = curve.type+"("+repr(curve.eid)+") = {"
    for p in range(len(curve.points)-1):
      line += repr(curve.points[p].eid)+", "
    line += repr(curve.points[-1].eid)+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += os.linesep
    self.lines.append(line)

  def addinterpcurve(self, interpcurve, comment=None):
    for curve in interpcurve.interpcurves: self.addcurve(curve, comment)

  def addsurface(self, surface, comment=None):
    if surface.eid:
      return
    surface.eid = self.sindex
    self.sindex += 1
    for curve in surface.curves:
      self.addcurve(curve)
    line = "Line Loop("+repr(self.cindex)+") = {"
    for c in range(len(surface.curves)-1):
      line += repr(surface.directions[c]*surface.curves[c].eid)+", "
    line += repr(surface.directions[-1]*surface.curves[-1].eid)+"};"+os.linesep
    self.lines.append(line)
    line = "Plane Surface("+repr(surface.eid)+") = {"+repr(self.cindex)+"};"
    if comment:
      line += " "+self.comment(comment)
    else:
      line += os.linesep
    self.lines.append(line)
    self.cindex += 1

  def addphysicalpoint(self, pid, points):
    line = "Physical Point("+repr(pid)+") = {"
    for p in range(len(points)-1):
      line += repr(points[p].eid)+", "
    line += repr(points[-1].eid)+"};"+os.linesep
    self.lines.append(line)

  def addphysicalline(self, pid, curves):
    line = "Physical Line("+repr(pid)+") = {"
    for c in range(len(curves)-1):
      line += repr(curves[c].eid)+", "
    line += repr(curves[-1].eid)+"};"+os.linesep
    self.lines.append(line)

  def addphysicalsurface(self, pid, surfaces):
    line = "Physical Surface("+repr(pid)+") = {"
    for s in range(len(surfaces)-1):
      line += repr(surfaces[s].eid)+", "
    line += repr(surfaces[-1].eid)+"};"+os.linesep
    self.lines.append(line)

  def addembed(self, surface, items):
    line = items[0].type+" {"
    for i in range(len(items)-1):
      assert(items[i].type==items[0].type)
      line += repr(items[i].eid)+", "
    line += repr(items[-1].eid)+"} In Surface {"+repr(surface.eid)+"};"+os.linesep
    self.lines.append(line) 

  def addtransfinitecurve(self, curves, n):
    for c in range(len(curves)):
      if curves[c].eid is None:
        self.addcurve(curves[c])
    line = "Transfinite Line {"
    for c in range(len(curves)-1):
      line += repr(curves[c].eid)+", "
    line += repr(curves[-1].eid)+"} = "+repr(n)+";"+os.linesep
    self.lines.append(line)

  def addtransfinitesurface(self, surface, corners, direction):
    if surface.eid is None:
      self.addsurface(surface)
    line = "Transfinite Surface {"+repr(surface.eid)+"} = {"
    for c in range(len(corners)-1):
      line += repr(corners[c].eid)+", "
    line += repr(corners[-1].eid)+"} "+direction+";"+os.linesep
    self.lines.append(line)

  def linebreak(self):
    self.lines.append(os.linesep)

  def comment(self, string):
    return "// "+string+os.linesep

  def produced_comment(self):
    self.lines.append(self.comment("Produced by: "+" ".join(sys.argv)))

  def write(self, filename):
    f = open(filename, 'w')
    f.writelines(self.lines)
    f.close()

class ElementaryEntity:
  def __init__(self, name=None):
    self.name = name
    self.eid  = None

  def cleareid(self):
    self.eid = None

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

  def cleareid(self):
    for point in self.points: point.cleareid()
    ElementaryEntity.cleareid(self)

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
    polys.append(interp.PiecewisePolynomial(u, list(zip(self.x, derivs[0])), orders=3, direction=None))
    polys.append(interp.PiecewisePolynomial(u, list(zip(self.y, derivs[1])), orders=3, direction=None))
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
      out = numpy.where(self.x < left)[0]
      if (len(out)>0):
        coord = self.intersectx(left)
        self.croppoint(out, coord)
    if right:
      out = numpy.where(self.x > right)[0]
      if (len(out)>0):
        coord = self.intersectx(right)
        self.croppoint(out, coord)
    if bottom:
      out = numpy.where(self.y < bottom)[0]
      if (len(out)>0):
        coord = self.intersecty(bottom)
        self.croppoint(out, coord)
    if top:
      out = numpy.where(self.y > top)[0]
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
    ind = numpy.where(self.x==xint)[0]
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpointy(self, yint):
    ind = numpy.where(self.y==yint)[0]
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
    pind = list(range(len(points)))
    derivs = self.find_derivatives()
    for p in ps:
      if p not in points or p==points[0] or p==points[-1]: continue
      for i in range(1,len(points)-1):
        if points[i] == p: break
      newp0 = Point([p.x-derivs[0][pind[i]], p.y-derivs[1][pind[i]]], res=p.res)
      newp1 = Point([p.x+derivs[0][pind[i]], p.y+derivs[1][pind[i]]], res=p.res)
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
        name = self.name+"_split"+repr(p)
      else: 
        name = None
      splines.append(Spline(newpoints, name=name, pid=self.pid))
    return splines

class InterpolatedCubicSpline:
  def __init__(self, points, name=None, pids=None, bctype='natural'):
    self.type = "InterpolatedCubicSpline"
    self.bctype=bctype
    self.points = [point for point in points]
    self.controlpoints = [point for point in points]
    self.controlpoints.sort(key=lambda point: point.x)
    self.name = name
    if pids: 
      assert(len(self.pids)==len(self.points)-1)
      self.pids = pids
    else:
      self.pids = [None for i in range(len(self.points)-1)]
    self.x = None
    self.y = None
    self.u = None
    self.interpu = None
    self.interpcurves = None
    self.length = None
    self.update()

  def __call__(self, delu, x0=None, der=0):
    # NOTE: this returns [x, dy/dx] when der=1...
    x = self.delu2x(delu, x0=x0)
    return [x, float(self.cs(x, nu=der))]

  def x2delu(self, x, x0=None):
    """Convert from Delta x to Delta u:
       u = \int_{x0}^{x} sqrt(1 + (dy(x')/dx')**2) dx'
       x0 is the lower bound of integration - provide to get an incremental u"""
    if x0 is None: x0 = self.x[0]
    return integ.quad(lambda xp: self.du(xp), x0, x)[0]/self.length

  def delu2x(self, delu, x0=None):
    """Convert from u to x:
       u(x) = \int_{x0}^x sqrt(1 + (dy(x')/dx')**2) dx'
       x0 is the lower bound of integration."""
    if x0 is None: x0 = self.x[0]
    return opt.fsolve(lambda x: self.x2delu(x, x0=x0)-delu, x0, fprime=lambda x: [self.du(x)])[0]

  def du(self, x):
    return sqrt(1.+float(self.cs(x, nu=1))**2)

  def update(self):
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    isort = numpy.argsort(self.x)
    points = numpy.asarray(self.points)[isort]
    self.points = points.tolist()
    pids = numpy.asarray(self.pids)[[i for i in isort if i != (len(isort)-1)]]
    self.pids = pids.tolist()
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    self.y = numpy.array([self.points[i].y for i in range(len(self.points))])
    controlx = numpy.array([cp.x for cp in self.controlpoints])
    controly = numpy.array([cp.y for cp in self.controlpoints])
    self.cs = interp.CubicSpline(controlx, controly, bc_type=self.bctype)
    self.length = 1.0 # must set this to one first to get the next line correct
    self.length = self.x2delu(self.x[-1])
    u = numpy.asarray([0.0])
    u = numpy.append(u, [self.x2delu(self.x[i], x0=self.x[i-1]) for i in range(1,len(self.x))])
    self.u = numpy.cumsum(u)
  
  def updateinterp(self):
    self.interpu = []
    self.interpcurves = []
    for p in range(len(self.points)-1):
      pid = self.pids[p]
      lengthfraction = self.u[p+1]-self.u[p]
      res0 = self.points[p].res/self.length/lengthfraction
      res1 = self.points[p+1].res/self.length/lengthfraction
      t = 0.0
      ts = [t]
      while t < 1.0:
        t = ts[-1] + (1.0 - t)*res0 + t*res1
        ts.append(t)
      ts = numpy.array(ts)/ts[-1]
      ls = (ts[1:]-ts[:-1])*lengthfraction*self.length
      res = [max(ls[i], ls[i+1]) for i in range(len(ls)-1)]
      point = self.points[p]
      self.interpu.append(self.u[p])
      for i in range(1,len(ts)-1):
        t = ts[i]
        self.interpu.append(self.u[p] + t*lengthfraction)
        npoint = Point(self(t*lengthfraction, x0=self.points[p].x), res=res[i-1])
        self.interpcurves.append(Line([point, npoint], pid))
        point = npoint
      npoint = self.points[p+1]
      self.interpcurves.append(Line([point, npoint], pid))
    self.interpu.append(self.u[p+1])
    self.interpu = numpy.asarray(self.interpu)

  def copyinterp(self, spline):
    assert(len(self.points)==len(spline.points))
    self.interpu = []
    self.interpcurves = []
    for p in range(len(self.points)-1):
      pid = self.pids[p]
      lengthfraction = self.u[p+1]-self.u[p]
      splinelengthfraction = spline.u[p+1]-spline.u[p]
      lengthratio = lengthfraction/splinelengthfraction
      splinepoint0 = spline.points[p]
      splinepoint1 = spline.points[p+1]
      splineus = spline.interpusinterval(splinepoint0, splinepoint1)
      splinecurves = spline.interpcurvesinterval(splinepoint0, splinepoint1)
      assert(len(splineus)==len(splinecurves)+1)
      point = self.points[p]
      self.interpu.append(self.u[p])
      for i in range(1, len(splineus)-1):
        self.interpu.append(self.u[p] + (splineus[i]-spline.u[p])*lengthratio)
        npoint = Point(self((splineus[i]-spline.u[p])*lengthratio, x0=self.points[p].x), \
                       res=splinecurves[i-1].points[1].res*lengthratio)
        self.interpcurves.append(Line([point,npoint], pid))
        point = npoint
      npoint = self.points[p+1]
      self.interpcurves.append(Line([point, npoint], pid))
    self.interpu.append(self.u[p+1])
    self.interpu = numpy.asarray(self.interpu)

  def intersecty(self, yint, extrapolate=False):
    return [self.cs.solve(yint, extrapolate=extrapolate)[0], yint]
      
  def intersectx(self, xint):
    return [xint, float(self.cs(xint))]

  def interpcurveindex(self, u):
    loc = abs(self.interpu - u).argmin()
    if loc == 0: 
      loc0 = loc
    elif loc == len(self.interpu)-1: 
      loc0 = loc-1
    else:
      if self.interpu[loc] < u: 
        loc0 = loc
      else: 
        loc0 = loc-1
    return loc0

  def unittangentx(self, x):
    der = self.cs(x, nu=1)
    vec = numpy.array([1.0, der])
    vec = vec/sqrt(sum(vec**2))
    return vec

  def translatenormal(self, dist):
    for p in range(len(self.points)):
      der = float(self.cs(self.points[p].x, nu=1))
      vec = numpy.array([-der, 1.0])
      vec = vec/sqrt(sum(vec**2))
      self.points[p].x += dist*vec[0]
      self.points[p].y += dist*vec[1]
    self.update()

  def croppoint(self, pind, coord):
    i = -1
    for i in range(len(pind)-1):
      p = self.points.pop(pind[i])
      if p in self.controlpoints:
        self.controlpoints.pop(self.controlpoints.index(p))
    self.points[pind[-1] - (i+1)].x = coord[0]
    self.points[pind[-1] - (i+1)].y = coord[1]
    self.update()

  def crop(self, left=None, bottom=None, right=None, top=None):
    if left is not None:
      out = numpy.where(self.x < left)[0]
      if (len(out)>0):
        coord = self.intersectx(left)
        self.croppoint(out, coord)
    if right is not None:
      out = numpy.where(self.x > right)[0]
      if (len(out)>0):
        coord = self.intersectx(right)
        self.croppoint(out, coord)
    if bottom is not None:
      out = numpy.where(self.y < bottom)[0]
      if (len(out)>0):
        coord = self.intersecty(bottom)
        self.croppoint(out, coord)
    if top is not None:
      out = numpy.where(self.y > top)[0]
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
    ind = numpy.where(self.x==xint)[0]
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpointy(self, yint):
    ind = numpy.where(self.y==yint)[0]
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpoint(self, name):
    for p in self.points:
      if p.name == name: return p
    return None

  def findpointindex(self, name):
    for p in range(len(self.points)):
      if self.points[p].name == name: return p
    return None

  def interpcurvesinterval(self, point0, point1):
    for l0 in range(len(self.interpcurves)):
      if self.interpcurves[l0].points[0] == point0: break

    for l1 in range(l0, len(self.interpcurves)):
      if self.interpcurves[l1].points[1] == point1: break

    return self.interpcurves[l0:l1+1]

  def interpusinterval(self, point0, point1):
    for l0 in range(len(self.interpcurves)):
      if self.interpcurves[l0].points[0] == point0: break

    for l1 in range(l0, len(self.interpcurves)):
      if self.interpcurves[l1].points[1] == point1: break

    return self.interpu[l0:l1+2]

  def appendpoint(self, p, pid=None):
    self.points.append(p)
    self.pids.append(pid)
    self.update()

  def cleareid(self):
    for curve in self.interpcurves: curve.cleareid()

class InterpolatedSciPySpline:
  def __init__(self, points, name=None, pids=None):
    self.type = "InterpolatedSciPySpline"
    self.points = points
    self.name = name
    if pids: 
      assert(len(self.pids)==len(self.points)-1)
      self.pids = pids
    else:
      self.pids = [None for i in range(len(self.points)-1)]
    self.x = None
    self.y = None
    self.u = None
    self.tck = None
    self.interpu = None
    self.interpcurves = None
    self.length = None
    self.update()

  def __call__(self, u, der=0):
    x, y = interp.splev(u, self.tck, der=der)
    return [float(x), float(y)]

  def update(self):
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    isort = numpy.argsort(self.x)
    points = numpy.asarray(self.points)[isort]
    self.points = points.tolist()
    pids = numpy.asarray(self.pids)[[i for i in isort if i != (len(isort)-1)]]
    self.pids = pids.tolist()
    self.x = numpy.array([self.points[i].x for i in range(len(self.points))])
    self.y = numpy.array([self.points[i].y for i in range(len(self.points))])
    self.tck, self.u = interp.splprep([self.x,self.y], s=0)
    self.length = integ.quad(lambda x: sqrt(sum(numpy.array(self(x, der=1))**2)), 0., 1.)[0]
  
  def updateinterp(self):
    self.interpu = []
    self.interpcurves = []
    for p in range(len(self.points)-1):
      pid = self.pids[p]
      lengthfraction = self.u[p+1]-self.u[p]
      res0 = self.points[p].res/self.length/lengthfraction
      res1 = self.points[p+1].res/self.length/lengthfraction
      t = 0.0
      ts = [t]
      while t < 1.0:
        t = ts[-1] + (1.0 - t)*res0 + t*res1
        ts.append(t)
      ts = numpy.array(ts)/ts[-1]
      ls = (ts[1:]-ts[:-1])*lengthfraction*self.length
      res = [max(ls[i], ls[i+1]) for i in range(len(ls)-1)]
      point = self.points[p]
      self.interpu.append(self.u[p])
      for i in range(1,len(ts)-1):
        t = ts[i]
        self.interpu.append(self.u[p] + t*lengthfraction)
        npoint = Point(self(self.interpu[-1]), res=res[i-1])
        self.interpcurves.append(Line([point, npoint], pid))
        point = npoint
      npoint = self.points[p+1]
      self.interpcurves.append(Line([point, npoint], pid))
    self.interpu.append(self.u[p+1])
    self.interpu = numpy.asarray(self.interpu)

  def copyinterp(self, spline):
    assert(len(self.points)==len(spline.points))
    self.interpu = []
    self.interpcurves = []
    for p in range(len(self.points)-1):
      pid = self.pids[p]
      lengthfraction = self.u[p+1]-self.u[p]
      splinelengthfraction = spline.u[p+1]-spline.u[p]
      lengthratio = lengthfraction/splinelengthfraction
      splinepoint0 = spline.points[p]
      splinepoint1 = spline.points[p+1]
      splineus = spline.interpusinterval(splinepoint0, splinepoint1)
      splinecurves = spline.interpcurvesinterval(splinepoint0, splinepoint1)
      assert(len(splineus)==len(splinecurves)+1)
      point = self.points[p]
      self.interpu.append(self.u[p])
      for i in range(1, len(splineus)-1):
        self.interpu.append(self.u[p] + (splineus[i]-spline.u[p])*lengthratio)
        npoint = Point(self(self.interpu[-1]), res=splinecurves[i-1].points[1].res*lengthratio)
        self.interpcurves.append(Line([point,npoint], pid))
        point = npoint
      npoint = self.points[p+1]
      self.interpcurves.append(Line([point, npoint], pid))
    self.interpu.append(self.u[p+1])
    self.interpu = numpy.asarray(self.interpu)

  def copytck(self):
    tck = []
    tck.append(self.tck[0].copy())
    tck.append([self.tck[1][0].copy(), self.tck[1][1].copy()])
    tck.append(self.tck[2])
    return tck

  def uintersectxy(self, xint, yint):
    tck = self.copytck()
    tck[1][0] = tck[1][0]-xint
    tck[1][1] = tck[1][1]-yint

    if tck[1][0][0]==0.0:  return [[0.0], []]
    if tck[1][0][-1]==0.0: return [[1.0], []]
    if tck[1][1][0]==0.0:  return [[], [0.0]]
    if tck[1][1][-1]==0.0: return [[], [1.0]]

    if numpy.all(numpy.sign(tck[1][0])==numpy.sign(tck[1][0])[0]) \
       and numpy.all(numpy.sign(tck[1][1])==numpy.sign(tck[1][1])[0]):
      if numpy.abs(tck[1][0]).min() < numpy.abs(tck[1][1]).min():
        assert(numpy.abs(tck[1][0]).min() < self.length*1.e-10)
        tck[1][0] = tck[1][0] - tck[1][0][numpy.abs(tck[1][0]).argmin()]
      else:
        assert(numpy.abs(tck[1][1]).min() < self.length*1.e-10)
        tck[1][1] = tck[1][1] - tck[1][1][numpy.abs(tck[1][1]).argmin()]

    t = tck[0]
    c0 = numpy.concatenate( (tck[1][0], numpy.array([0]*(len(t)-len(tck[1][0])), dtype=tck[1][0].dtype)))
    c1 = numpy.concatenate( (tck[1][1], numpy.array([0]*(len(t)-len(tck[1][1])), dtype=tck[1][1].dtype)))
    z0, m0,ierr = interp.fitpack2.dfitpack.sproot(t, c0)
    z1, m1,ierr = interp.fitpack2.dfitpack.sproot(t, c1)
    return [z0[:m0], z1[:m1]]
      
  def intersecty(self, yint):
    spoint = self.uintersecty(yint)
    assert(len(spoint)==1)
    x,y = interp.splev(spoint[0], self.tck)
    return [float(x),float(y)]
      
  def uintersecty(self, yint):
    tck = self.copytck()
    tck[1][1] = tck[1][1]-yint
    t = tck[0]
    c = numpy.concatenate( (tck[1][1], numpy.array([0]*(len(t)-len(tck[1][1])), dtype=tck[1][1].dtype)))
    z,m,ierr = interp.fitpack2.dfitpack.sproot(t, c)
    return z[:m]
      
  def intersectx(self, xint):
    spoint = self.uintersectx(xint)
    assert(len(spoint)==1)
    x,y = interp.splev(spoint[0], self.tck)
    return [float(x),float(y)]

  def uintersectx(self, xint):
    tck = self.copytck()
    tck[1][0] = tck[1][0]-xint
    t = tck[0]
    c = numpy.concatenate( (tck[1][0], numpy.array([0]*(len(t)-len(tck[1][0])), dtype=tck[1][0].dtype)))
    z,m,ierr = interp.fitpack2.dfitpack.sproot(t, c)
    return z[:m]

  def uxymindist(self, x, y):
    spoint = self.uintersectxy(x,y)
    for i in range(2): 
      if len(spoint[i])>0: break
    f = lambda u: numpy.sqrt(sum((numpy.asarray(self(u)) - numpy.asarray([x,y]))**2))
    return opt.fmin_bfgs(f, spoint[i][0], disp=0)[0]

  def uxyapprox(self, x, y):
    spoint = self.uintersectxy(x, y)
    for i in range(2): 
      if len(spoint[i])>0: break
    if len(spoint[i])==0: print("x, y = ", x, y)
    loc = self.interpcurveindex(spoint[i][0])
    totallength = sqrt((self.interpcurves[loc].points[-1].x - self.interpcurves[loc].points[0].x)**2 \
                      +(self.interpcurves[loc].points[-1].y - self.interpcurves[loc].points[0].y)**2)
    partiallength = sqrt((x - self.interpcurves[loc].points[0].x)**2 \
                        +(y - self.interpcurves[loc].points[0].y)**2)
    fractionallength = partiallength/totallength
    return self.interpu[loc] + fractionallength*(self.interpu[loc+1]-self.interpu[loc])

  def interpcurveindex(self, u):
    loc = abs(self.interpu - u).argmin()
    if loc == 0: 
      loc0 = loc
    elif loc == len(self.interpu)-1: 
      loc0 = loc-1
    else:
      if self.interpu[loc] < u: 
        loc0 = loc
      else: 
        loc0 = loc-1
    return loc0

  def unittangent(self, u):
    der = self(u, der=1)
    vec = numpy.array([der[0], der[1]])
    vec = vec/sqrt(sum(vec**2))
    return vec

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
    if left is not None:
      out = numpy.where(self.x < left)[0]
      if (len(out)>0):
        coord = self.intersectx(left)
        self.croppoint(out, coord)
    if right is not None:
      out = numpy.where(self.x > right)[0]
      if (len(out)>0):
        coord = self.intersectx(right)
        self.croppoint(out, coord)
    if bottom is not None:
      out = numpy.where(self.y < bottom)[0]
      if (len(out)>0):
        coord = self.intersecty(bottom)
        self.croppoint(out, coord)
    if top is not None:
      out = numpy.where(self.y > top)[0]
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
    ind = numpy.where(self.x==xint)[0]
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpointy(self, yint):
    ind = numpy.where(self.y==yint)[0]
    if (len(ind)==0): 
      return None
    else:
      return self.points[ind[0]]
  
  def findpoint(self, name):
    for p in self.points:
      if p.name == name: return p
    return None

  def findpointindex(self, name):
    for p in range(len(self.points)):
      if self.points[p].name == name: return p
    return None

  def interpcurvesinterval(self, point0, point1):
    for l0 in range(len(self.interpcurves)):
      if self.interpcurves[l0].points[0] == point0: break

    for l1 in range(l0, len(self.interpcurves)):
      if self.interpcurves[l1].points[1] == point1: break

    return self.interpcurves[l0:l1+1]

  def interpusinterval(self, point0, point1):
    for l0 in range(len(self.interpcurves)):
      if self.interpcurves[l0].points[0] == point0: break

    for l1 in range(l0, len(self.interpcurves)):
      if self.interpcurves[l1].points[1] == point1: break

    return self.interpu[l0:l1+2]

  def appendpoint(self, p, pid=None):
    self.points.append(p)
    self.pids.append(pid)
    self.update()

  def cleareid(self):
    for curve in self.interpcurves: curve.cleareid()

class Surface(ElementaryEntity, PhysicalEntity):
  def __init__(self, curves, name=None, pid = None):
    self.type = "Surface"
    self.curves = [curves[0]]
    self.directions = [1]
    ind = list(range(1, len(curves)))
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
      if i == len(ind)-1 and len(self.directions) == j:
        print("Failed to complete line loop.")
        sys.exit(1)
      self.curves.append(curves[ind[i]])
      i = ind.pop(i)
    ElementaryEntity.__init__(self, name=name)
    PhysicalEntity.__init__(self, pid=pid)

  def cleareid(self):
    for curve in self.curves: curve.cleareid()
    ElementaryEntity.cleareid(self)

class Geometry:
  def __init__(self):
    self.namedpoints = {}
    self.namedcurves = {}
    self.namedinterpcurves = {}
    self.namedsurfaces = {}
    self.physicalpoints = {}
    self.physicalcurves = {}
    self.physicalsurfaces = {}
    self.points  = []
    self.curves  = []
    self.interpcurves  = []
    self.surfaces = []
    self.pointcomments  = []
    self.curvecomments  = []
    self.interpcurvecomments  = []
    self.surfacecomments = []
    self.pointembeds = {}
    self.lineembeds = {}
    self.transfinitecurves = {}
    self.transfinitesurfaces = {}
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

  def addtransfinitecurve(self, curve, name=None, comment=None, n=2):
    self.addcurve(curve, name, comment)
    if n in self.transfinitecurves:
      self.transfinitecurves[n].append(curve)
    else:
      self.transfinitecurves[n] = [curve]

  def addtransfinitesurface(self, surface, corners, direction="Right", name=None, comment=None):
    self.addsurface(surface, name, comment)
    self.transfinitesurfaces[surface] = (corners, direction)

  def addinterpcurve(self, interpcurve, name=None, comment=None):
    self.interpcurves.append(interpcurve)
    self.interpcurvecomments.append(comment)
    if not comment:
      if name:
        self.interpcurvecomments[-1] = name
      elif interpcurve.name:
        self.interpcurvecomments[-1] = interpcurve.name
    if name:
      self.namedinterpcurves[name] = interpcurve
    elif interpcurve.name:
      self.namedinterpcurves[interpcurve.name] = interpcurve
    for curve in interpcurve.interpcurves: self.addtransfinitecurve(curve)

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
    self.cleareid()
    
    self.geofile = GeoFile()
    for s in range(len(self.surfaces)):
      surface = self.surfaces[s]
      self.geofile.addsurface(surface, self.surfacecomments[s])
    self.geofile.linebreak()
    for c in range(len(self.curves)):
      curve = self.curves[c]
      self.geofile.addcurve(curve, self.curvecomments[c])
    self.geofile.linebreak()
    for p in range(len(self.points)):
      point = self.points[p]
      self.geofile.addpoint(point, self.pointcomments[p])
    self.geofile.linebreak()

    for surface,items in iter(sorted(list(self.pointembeds.items()), key=lambda item: item[0].eid)):
      self.geofile.addembed(surface,items)
    for surface,items in iter(sorted(list(self.lineembeds.items()), key=lambda item: item[0].eid)):
      self.geofile.addembed(surface,items)
    self.geofile.linebreak()

    for n,curves in iter(sorted(self.transfinitecurves.items())):
      self.geofile.addtransfinitecurve(curves,n)
    self.geofile.linebreak()
    for surface,corners in iter(sorted(list(self.transfinitesurfaces.items()), key=lambda item: item[0].eid)):
      self.geofile.addtransfinitesurface(surface,corners[0],corners[1])
    self.geofile.linebreak()
    
    for pid,surfaces in iter(sorted(self.physicalsurfaces.items())):
      self.geofile.addphysicalsurface(pid,surfaces)
    self.geofile.linebreak()
    for pid,curves in iter(sorted(self.physicalcurves.items())):
      self.geofile.addphysicalline(pid,curves)
    self.geofile.linebreak()
    for pid,points in iter(sorted(self.physicalpoints.items())):
      self.geofile.addphysicalpoint(pid,points)
    self.geofile.linebreak()

    self.geofile.produced_comment()

    self.geofile.write(filename)

  def pylabplot(self, lineres=100):
    import pylab
    for interpcurve in self.interpcurves:
      unew = numpy.arange(interpcurve.u[0], interpcurve.u[-1]+((interpcurve.u[-1]-interpcurve.u[0])/(2.*lineres)), 1./lineres)
      pylab.plot(interpcurve(unew)[0], interpcurve(unew)[1], 'b')
      pylab.plot(interpcurve.x, interpcurve.y, 'ob')
      for curve in interpcurve.interpcurves:
        unew = numpy.arange(curve.u[0], curve.u[-1]+((curve.u[-1]-curve.u[0])/(2.*lineres)), 1./lineres)
        pylab.plot(curve(unew)[0], curve(unew)[1], 'k')
        pylab.plot(curve.x, curve.y, '+k')
    for curve in self.curves:
      unew = numpy.arange(curve.u[0], curve.u[-1]+((curve.u[-1]-curve.u[0])/(2.*lineres)), 1./lineres)
      pylab.plot(curve(unew)[0], curve(unew)[1], 'k')
      pylab.plot(curve.x, curve.y, 'ok')
    #for point in self.points:
    #  pylab.plot(point.x, point.y, 'ok')

  def cleareid(self):
    for surface in self.surfaces: surface.cleareid()
    for interpcurve in self.interpcurves: interpcurve.cleareid()
    for curve in self.curves: curve.cleareid()
    for point in self.points: point.cleareid()
  
