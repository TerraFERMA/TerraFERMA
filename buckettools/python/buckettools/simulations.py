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

import os
import subprocess
import sys
import hashlib
import shutil
import threading
import Queue
import copy
import string
import glob
import re
import collections
from buckettools.threadlibspud import *
import traceback

####################################################################################

class SimulationsErrorInitialization(Exception):
  def __init__(self, msg = "Error while initializing simulation."):
    self.message = msg

class SimulationsErrorWriteOptions(Exception):
  def __init__(self, msg = "Error while writing options."):
    self.message = msg

class SimulationsErrorConfigure(Exception):
  def __init__(self, msg = "Error while configuring the simulation."):
    self.message = msg

class SimulationsErrorBuild(Exception):
  def __init__(self, msg = "Error while building the simulation."):
    self.message = msg

class SimulationsErrorRun(Exception):
  def __init__(self, msg = "Error while running simulation."):
    self.message = msg

class SimulationsErrorVariable(Exception):
  def __init__(self, msg = "Error while evaluating variables."):
    self.message = msg

class SimulationsErrorTest(Exception):
  def __init__(self, msg = "Error while testing simulations."):
    self.message = msg

####################################################################################

class ThreadIterator(list):
  '''A thread-safe iterator over a list.'''
  def __init__(self, seq):
    self.list=list(seq)
    self.lock=threading.Lock()

  def __iter__(self):
    return self

  def next(self):
    self.lock.acquire()
    if len(self.list)>0:
      ans=self.list.pop()
    else:
      ans=None
    self.lock.release()
    if ans is None:
      raise StopIteration
    return ans

####################################################################################

class TestOrVariable:
    """Tests and variables have a lot in common. This code unifies the commonalities."""
    def __init__(self, name, code):
        self.name = name
        self.code = code

####################################################################################

class Test(TestOrVariable):
    """A test for the model output"""
    def run(self, varsdict):
        tmpdict = copy.copy(varsdict) # don't let the test code modify the variables
        try:
          exec self.code in tmpdict
          return True
        except AssertionError:
          # in case of an AssertionError, we assume the test has just failed
          return False
        except:
          # tell us what else went wrong:
          traceback.print_exc()
          return False

####################################################################################

class Variable(TestOrVariable):
    """A variable definition for use in tests"""
    def run(self, varsdict):
        try:
            exec self.code in varsdict
        except:
            print "Variable computation raised an exception"
            print "-" * 80
            for (lineno, line) in enumerate(self.code.split('\n')):
              print "%3d  %s" % (lineno+1, line)
            print "-" * 80
            traceback.print_exc()
            print "-" * 80
            raise Exception

        if self.name not in varsdict.keys():
            print "self.name == ", self.name
            print "varsdict.keys() == ", varsdict.keys()
            print "self.name not found: does the variable define the right name?"
            raise Exception

####################################################################################

class NestedList(list):
  '''A class that implements a nested list structure based on a dictionary of paramters.  
     The resulting nested list can be indexed by a dictionary of parameters that works out which
     indices to index into the nested list structure.'''

  def __init__(self, param=None):
    '''Initialize a nested list.'''
    # initialize the base list class
    super(NestedList, self).__init__()

    # record the parameter dictionary or set up an empty one
    if param is None: param = collections.OrderedDict()
    self.param = param

    # if the parameters have been provided in an ordered dictionary then
    # we respect that order, otherwise we sort the keys alphabetically:
    if isinstance(self.param, collections.OrderedDict):
      self.sortedkeys = self.param.keys()
    else:
      self.sortedkeys = sorted(self.param.keys())

    # create the nested list structure, beneath the top level list object, self:
    self.createnestedlist(self)

  def createnestedlist(self, parentlist, level=0):
    '''Recursively set up a nested list.'''
    # if we've reached the bottom of the dictionary then append Nones and done recurse
    if level == len(self.sortedkeys)-1:
      for val in self.param[self.sortedkeys[level]]: parentlist.append(None)
    # if we haven't reached the bottom then loop over the values and recursively call
    # this function to set up the sublevels of the nested list
    elif level < len(self.sortedkeys)-1:
      for val in self.param[self.sortedkeys[level]]: 
        childlist = []
        self.createnestedlist(childlist, level=level+1)
        parentlist.append(childlist)

  def getindexdict(self, index):
    '''Based on an input index dictionary setup a full dictionary to index to nested list 
       that includes all keys and whose values are the indicies for the list (rather than 
       the string values we want).'''
    indexdict = {}
    # loop over all the keys we want
    for key in self.sortedkeys:
      # get all the values that are allowed for this key
      paramvals = self.param[key]
      # if the key is in the input then...
      if key in index:
        # it it has been mentioned then we only take the values requested
        indexvals = index[key] # which values are requested?
        if not isinstance(indexvals, list): indexvals = [indexvals] # in case we didn't ask for a list
        indicies = [] # list the indicies to the parameter dictionary
        for indexval in indexvals:
          try: 
            indicies.append(paramvals.index(indexval))
          except ValueError:
            # value we don't know about requested
            pass
      else:
        # if it hasn't been specifically mentioned then we take all values at that level
        indicies = range(len(paramvals))
      # for the current key record the indicies we want to fetch
      indexdict[key] = indicies
    # return the full indexdict
    return indexdict

  def slicenestedlist(self, parentlist, level=0):
    '''Recursively slice the nested list to find a particular value (in the indexdict) based on indexdict (which must be set).'''
    # if we've reached the bottom level then we can actually fetch the items
    if level == len(self.sortedkeys)-1:
      newlist = [parentlist[i] for i in self.indexdict[self.sortedkeys[level]]]
    # if we've not reached the bottom level then we need to recurse
    elif level < len(self.sortedkeys)-1:
      newlist = [self.slicenestedlist(parentlist[i], level=level+1) for i in self.indexdict[self.sortedkeys[level]]]

    # collapse lists of length 1 so we don't end up with unnecessary levels in the returned item
    # NOTE: good idea or not?
    if len(newlist)==1:
      return newlist[0]
    else:
      return newlist

  def __getitem__(self, index):
    '''Fetch item from nested list.  index can be dictionary that references parameters dictionary or standard integer.'''
    # if we have a dictionary index
    if isinstance(index, dict):
      self.indexdict = self.getindexdict(index) # set up the full indexdict (necessary before next call)
      item = self.slicenestedlist(self) # fetch the item(s) by slicing the nested list
      del self.indexdict # delete it again so that we don't get confused on the next call
      return item
    # else assume this is a standard integer index and call the base list object
    else:
      return super(NestedList, self).__getitem__(index)

  def setnestedlistvalue(self, parentlist, value, level=0):
    '''Recursively set a value in the nested list based on the indexing dictionary.'''
    # if we've reached the bottom then set the value
    if level == len(self.sortedkeys)-1:
      for i in self.indexdict[self.sortedkeys[level]]: parentlist[i] = value
    # otherwise recurse to the next level
    elif level < len(self.sortedkeys)-1:
      for i in self.indexdict[self.sortedkeys[level]]: self.setnestedlistvalue(parentlist[i], value, level=level+1)

  def __setitem__(self, index, value):
    '''Set item in nested list.  index can be dictionary that references parameters dictionary or standard integer.'''
    # if we have a ductionary index
    if isinstance(index, dict):
      self.indexdict = self.getindexdict(index) # set up the full indexdict (necessary before next call)
      self.setnestedlistvalue(self, value) # set the item(s) by recursing through the nested list
      del self.indexdict # delete the indexdict so that we don't get confused on the next call
      return None
    # else assume this is a standard integer index and call the base list object
    else:
      return super(NestedList, self).__setitem__(index, value)

####################################################################################

class Run:
  def __init__(self, path, optionsdict, currentdirectory, \
               logprefix=None, dependents=[]):
    self.optionsdict = optionsdict

    filename = os.path.basename(path)
    self.filename, self.ext = os.path.splitext(filename)

    self.name = self.optionsdict["name"]

    self.logprefix = logprefix

    self.currentdirectory = currentdirectory

    # we preserve the original path passed to us so we can use it to reference back to this
    # simulation in the batch global options dictionary
    self.path = path
    # but we also want to make sure that we have the correct base directory to work from
    dirname = os.path.dirname(path)
    if os.path.isabs(dirname):
      self.basedirectory = os.path.normpath(dirname)
    else:
      self.basedirectory = os.path.normpath(os.path.join(self.currentdirectory, dirname))
    self.rundirectory = self.getrundirectory()
    self.runinputdirectory = self.basedirectory

    if "nruns" in self.optionsdict:
      self.nruns = self.optionsdict["nruns"]
    else:
      self.nruns = 1
    self.dependents = dependents
    self.lock=threading.Lock()

    self.variables = [Variable(name, code) for name, code in self.optionsdict["variables"].iteritems()]

    self.dependencies = []
    self.alreadyrun = False

    self.deferred = False # can't have deferred runs when you have no dependencies

  def writeoptions(self):
    
    inputfile = open(os.path.join(self.runinputdirectory, self.filename+self.ext))
    inputstr = inputfile.read()
    inputfile.close()
    
    valuesdict = {"input_file":inputstr}
    self.updateoptions(valuesdict=valuesdict)

    try:
      os.makedirs(self.rundirectory)
    except OSError:
      pass
    
    inputstr = valuesdict["input_file"]
    outputfile = open(os.path.join(self.rundirectory, self.filename+self.ext), 'w')
    outputfile.write(inputstr)
    outputfile.close()

  def updateoptions(self, valuesdict=None, prefix=""):
    os.chdir(self.basedirectory)
    if valuesdict is None: valuesdict={}
    for k,v in self.optionsdict[prefix+"values"].iteritems():
      if k in valuesdict:
        self.log("ERROR: in updateoptions, %s multiply defined"%(k))
        if k=="input_file":
          self.log("parameter cannot be named input_file")
        raise SimulationsErrorWriteOptions
      valuesdict[k] = v
    for update in self.optionsdict[prefix+"updates"].itervalues():
      if update is not None:
        exec update in valuesdict
    os.chdir(self.currentdirectory)

  def configure(self, force=False):
    return None

  def build(self, force=False):
    return None

  def run(self, force=False):
    
    self.lock.acquire()

    error = False
    
    for dependency in self.dependencies: dependency.run(force=force)

    self.log("checking in directory %s"%(os.path.relpath(self.rundirectory, self.currentdirectory)))
    if not self.alreadyrun:
      commands = self.getcommands()

      for r in xrange(self.nruns):
        requiredinput = self.getrequiredinput(r)
        requiredoutput = self.getrequiredoutput(r)

        dirname = os.path.join(self.rundirectory, "run_"+`r`.zfill(len(`self.nruns`)))
        try:
          os.makedirs(dirname)
        except OSError:
          pass

        input_changed = False
        for filepath in requiredinput:
          try:
            checksum = hashlib.md5(open(os.path.join(dirname, os.path.basename(filepath))).read()).hexdigest()
          except:
            checksum = None
          try:
            input_changed = input_changed or checksum != hashlib.md5(open(filepath).read()).hexdigest()
          except IOError:
            self.log("WARNING: Unable to open %s"%(filepath))
            input_changed = True

        output_missing = False
        for filepath in requiredoutput:
          try:
            output_file = open(os.path.join(dirname, filepath))
            output_file.close()
          except IOError:
            output_missing = True
        if len(requiredoutput)==0: output_missing=True # don't know what output is needed so we have to force
          
        if output_missing or input_changed or force:
          self.log("  running in directory %s"%(os.path.relpath(dirname, self.currentdirectory)))
          # file has changed or a recompilation is necessary
          for filepath in requiredoutput:
            try:
              os.remove(os.path.join(dirname, filepath))
            except OSError:
              pass
          for filepath in requiredinput:
            try:
              shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
            except IOError:
              self.log("WARNING: required input (%s) not found, continuing anyway."%(filepath))
          
          env = copy.copy(os.environ)
          try:
            env["PYTHONPATH"] = ":".join([dirname, env["PYTHONPATH"]])
          except KeyError:
            env["PYTHONPATH"] = dirname

          for command in commands:
            p = subprocess.Popen(command, cwd=dirname, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
            retcode = p.wait()
            if retcode!=0:
              self.log("ERROR: Command %s returned %d in directory: %s"%(" ".join(command), retcode, dirname))
              error = True
              break
          
          self.log("  finished in directory %s:"%(os.path.relpath(dirname, self.currentdirectory)))

      self.alreadyrun = True
    self.lock.release()
    if error: raise SimulationsErrorRun

  def clean(self):
    try:
      shutil.rmtree(os.path.join(self.basedirectory, self.filename+self.ext+".run"))
    except OSError:
      pass

  def evaluatevariables(self):
    '''Evaluate the variables for this simulation and its dependencies.'''

    error = False
    # set up the variables dictionary
    varsdict = {}
    # loop over all the runs
    for r in xrange(self.nruns):
      # change to the run directory to run the python variable assignment
      dirname = os.path.join(self.rundirectory, "run_"+`r`.zfill(len(`self.nruns`)))
      try:
        os.chdir(dirname)
      except OSError:
        self.log("WARNING: could not change to directory %s"%(dirname))

      # loop over the variables
      for var in self.variables:
        if self.nruns > 1: varsdict = []
        # set up a temporary variable dictionary for this variable calculation
        tmpdict  = {}
        # try running the variable assignment
        try:
          var.run(tmpdict)
        except:
          self.log("ERROR: failure while calculating variable %s." % (str(var.name)))
          error = True
          continue

        if self.nruns>1:
          varsdict[var.name].append(tmpdict[var.name])
        else:
          varsdict[var.name] = tmpdict[var.name]

        self.log("  Assigned %s = %s" % (str(var.name), str(varsdict[var.name])))
    # make sure we change back to the current directory
    os.chdir(self.currentdirectory)

    if error: raise SimulationsErrorVariable

    return varsdict

  def log(self, string):
    if self.logprefix is not None: string = self.logprefix+": "+string
    print string

  def getrundirectory(self):
    '''Return the run directory for this simulation.'''

    # start at the base directory
    rundirectory = os.path.join(self.basedirectory, self.filename+self.ext+".run")

    # preserve the order of an ordered dictionary or use alphabetical
    if isinstance(self.optionsdict["values"], collections.OrderedDict):
      sorteditems = self.optionsdict["values"].items()
    else:
      sorteditems = sorted(self.optionsdict["values"].items())

    # append directory levels for all parameter values
    for directory in [item[0]+"_"+item[1] for item in sorteditems]: 
      rundirectory = os.path.join(rundirectory, directory)

    # return the run directory
    return rundirectory

  def getdependencyoptions(self, options):
    suboptionsdict = collections.OrderedDict()
    for option in options: suboptionsdict[option] = self.optionsdict["values"][option]
    return suboptionsdict

  def getrequiredoutput(self, run):
    requiredoutput = []
    if "output" in self.optionsdict: requiredoutput += self.optionsdict["output"]
    return requiredoutput
 
  def getdependencyrequiredoutput(self, run=0):
    requiredoutput = []
    for dependency in self.dependencies:
      if self.nruns == dependency.nruns: 
        rundir = "run_"+`run`.zfill(len(`self.nruns`))
      else:
        rundir = "run_"+`0`.zfill(len(`dependency.nruns`))
      requiredoutput += [os.path.join(dependency.rundirectory, rundir, output) \
                                           for output in dependency.getrequiredoutput(run)]
      requiredoutput += dependency.getdependencyrequiredoutput(run)
    return requiredoutput

  def listinput(self):
    requiredinput = []
    if not self.deferred: requiredinput.append(os.path.join(self.basedirectory, self.filename+self.ext))
    if "input" in self.optionsdict:
      requiredinput += self.optionsdict["input"]
    return requiredinput

  def getrequiredinput(self, run):
    requiredinput = [os.path.join(self.rundirectory, self.filename+self.ext)]
    # get any input specified in the optionsdict
    if "input" in self.optionsdict: 
      for filepath in self.optionsdict["input"]:
        # filter out any filenames that have already been added (i.e. make sure our latest input file is used)
        if os.path.basename(filepath) not in [os.path.basename(inputpath) for inputpath in requiredinput]:
          requiredinput.append(filepath)
    # get any output from dependencies
    for filepath in self.getdependencyrequiredoutput(run):
      # filter out any filenames that have already been added (highest dependencies take priority over lower ones)
      if os.path.basename(filepath) not in [os.path.basename(inputpath) for inputpath in requiredinput]:
        requiredinput.append(filepath)
    return requiredinput

  def getcommands(self):
    commands = []
    if "run" in self.optionsdict: commands = self.optionsdict["run"]
    return commands

################################################################################################

class Simulation(Run):
  def __init__(self, path, optionsdict, currentdirectory, tfdirectory, \
               logprefix=None, dependents=[]):

    Run.__init__(self, path, optionsdict, currentdirectory, \
                 logprefix=logprefix, dependents=dependents)

    self.tfdirectory = tfdirectory
    self.builddirectory = self.getbuilddirectory()

    self.dependencies = []
    if "dependencies" in self.optionsdict:
       for dependencypath, depoptionsdict in self.optionsdict["dependencies"].iteritems():
         depoptionsdict["values"] = self.getdependencyoptions(depoptionsdict["updates"].keys())
         if depoptionsdict["type"] is Run:
           self.dependencies.append(Run(dependencypath, depoptionsdict, \
                                        currentdirectory, logprefix=logprefix, \
                                        dependents=[self]))
         elif depoptionsdict["type"] is Simulation:
           self.dependencies.append(Simulation(dependencypath, depoptionsdict, \
                                               currentdirectory, tfdirectory, \
                                               logprefix=logprefix, dependents=[self]))
         else:
           self.log("ERROR: Unknown dependency simulation type.")
           raise SimulationsErrorInitialization

    # if the input options file is an the output of a dependency then we have to
    # defer generating, configuring and building until runtime
    filepaths = [os.path.dirname(depoutput)  \
                 for depoutput in self.getdependencyrequiredoutput() \
                 if self.filename+self.ext == os.path.basename(depoutput)]
    self.deferred = len(filepaths)>0
    if self.deferred:
      if os.path.isabs(filepaths[0]):
        self.runinputdirectory = os.path.normpath(filepaths[0])
      else:
        self.runinputdirectory = os.path.normpath(os.path.join(currentdirectory, filepaths[0]))

  def writeoptions(self):
    
    threadlibspud.load_options(os.path.join(self.runinputdirectory, self.filename+self.ext))
    
    self.updateoptions()
    
    try:
      os.makedirs(self.rundirectory)
    except OSError:
      pass
    
    libspud.write_options(os.path.join(self.rundirectory, self.filename+self.ext))

    try:
      os.makedirs(self.builddirectory)
    except OSError:
      pass
    
    libspud.write_options(os.path.join(self.builddirectory, self.filename+self.ext))

    threadlibspud.clear_options()

  def writecheckpointoptions(self, basedir, basefile):
    threadlibspud.load_options(os.path.join(basedir, basefile))
    self.updateoptions(prefix="checkpoint_")
    try:
      os.makedirs(self.rundirectory)
    except OSError:
      pass
    libspud.write_options(os.path.join(basedir, "checkpoint", basefile))
    threadlibspud.clear_options()

  def configure(self, force=False):

    dirname = os.path.join(self.builddirectory, "build")

    if force: 
      try:
        shutil.rmtree(dirname)
      except OSError:
        pass
    
    try:
      os.makedirs(dirname)
    except OSError:
      pass

    p = subprocess.Popen(["cmake", \
                          "-DOPTIONSFILE="+os.path.join(self.builddirectory, self.filename+self.ext), \
                          "-DCMAKE_BUILD_TYPE=RelWithDebInfo", \
                          "-DLOGLEVEL=INFO", \
                          "-DEXECUTABLE="+self.filename, \
                          os.path.join(self.tfdirectory,os.curdir)],
                          cwd=dirname)
    retcode = p.wait()
    if retcode!=0:
      self.log("ERROR: cmake returned %d in directory: %s"%(retcode, dirname))
      raise SimulationsErrorConfigure

  def build(self, force=False):

    dirname = os.path.join(self.builddirectory, "build")
    if force:
      p = subprocess.Popen(["make", "clean"], cwd=dirname)
      retcode = p.wait()
      if retcode!=0:
        self.log("ERROR: make clean returned %d in directory: %s"%(retcode, dirname))
        raise SimulationsErrorBuild

    p = subprocess.Popen(["make"], cwd=dirname)
    retcode = p.wait()
    if retcode!=0:
      self.log("ERROR make returned %d in directory: %s"%(retcode, dirname))
      raise SimulationsErrorBuild

  def clean(self):
    Run.clean(self)

    try:
      shutil.rmtree(os.path.join(self.basedirectory, self.filename+self.ext+".build"))
    except OSError:
      pass

  def checkpointrun(self):

    if len(self.optionsdict["checkpoint_values"])==0: return
    
    for r in xrange(self.nruns):
      dirname = os.path.join(self.rundirectory, "run_"+`r`.zfill(len(`self.nruns`)))

      for root, dirnames, files in os.walk(dirname, topdown=True):
        depth = string.count(root, os.path.sep) - string.count(dirname, os.path.sep) + 1
        if "checkpoint" in dirnames:
          # there is a checkpoint directory beneath the current directory
          dirnames[:] = ["checkpoint"]
          files = glob.glob1(os.path.join(root, "checkpoint"), "*"+self.ext)
          files = [f for f in files if "_checkpoint"*(depth+1) in f]
          if len(files)==0: dirnames[:] = [] # but it contains no appropriate checkpoint files so don't descend
        else:
          # there is no checkpoint directory beneath the current directory
          dirnames[:] = []

      basedir = root

      self.log("checking in directory %s"%(os.path.relpath(basedir, self.currentdirectory)))

      files = glob.glob1(basedir, "*"+self.ext)
      files = [f for f in files if "_checkpoint"*depth in f]
      files = sorted(files, key=lambda f: int(f.split(self.ext)[0].split("_")[-1]))

      if len(files)==0: continue
      
      basefile = files[-1]

      requiredinput = self.getrequiredcheckpointinput(basedir, basefile)
      commands = self.getcheckpointcommands(basefile)

      dirname = os.path.join(basedir, "checkpoint")
      try:
        os.makedirs(dirname)
      except OSError:
        pass
        
      self.log("  running in directory %s"%(os.path.relpath(dirname, self.currentdirectory)))
      for filepath in requiredinput:
        shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
      self.writecheckpointoptions(basedir, basefile)
      
      for command in commands:
        p = subprocess.Popen(command, cwd=dirname)
        retcode = p.wait()
        if retcode!=0:
          self.log("ERROR Command %s returned %d in directory: %s"%(" ".join(command), retcode, dirname))
          raise SimulationsErrorRun
      
      self.log("  finished in directory %s"%(os.path.relpath(dirname, self.currentdirectory)))

  def getcheckpointcommands(self, basefile):
    commands = [[os.path.join(self.builddirectory, "build", self.filename), "-vINFO", "-l", basefile]]
    return commands

  def getcommands(self):
    commands = [[os.path.join(self.builddirectory, "build", self.filename), "-vINFO", "-l", self.filename+self.ext]]
    return commands

  def getbuilddirectory(self):
    '''Return the path to the build directory for this simulation.'''

    # start at the base directory
    builddirectory = os.path.join(self.basedirectory, self.filename+self.ext+".build")

    # decide the order (preserve what exists if ordered dictionary or use alphabetical)
    if isinstance(self.optionsdict["values"], collections.OrderedDict):
      sorteditems = self.optionsdict["values"].items()
    else:
      sorteditems = sorted(self.optionsdict["values"].items())

    # append directory levels if this is not a single build parameter
    for directory in [item[0]+"_"+item[1] for item in sorteditems if not self.optionsdict["builds"][item[0]]]:
      builddirectory = os.path.join(builddirectory, directory)

    # return the build directory
    return builddirectory

################################################################################################

class SimulationBatch:
  def __init__(self, globaloptionsdict, filename, currentdirectory, tfdirectory, \
               tests={}, nthreads=1):

    self.globaloptionsdict = globaloptionsdict
    self.currentdirectory = currentdirectory
    self.filename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    if os.path.isabs(dirname):
      self.basedirectory = os.path.normpath(dirname)
    else:
      self.basedirectory = os.path.normpath(os.path.join(currentdirectory, dirname))
    self.logprefix = os.path.relpath(filename, currentdirectory)
    self.nthreads = nthreads

    # set up the list of simulations/runs in this batch
    self.simulations = []
    # loop over all simulation paths and their global options dictionaries
    for simulationpath, simoptionsdict in self.globaloptionsdict.iteritems():
      # set up a list to contain the options dictionaries for each parameter set
      slicedoptionslist = []
      self.sliceoptionsdict(slicedoptionslist, simoptionsdict["values"]) # slice the global options dictionary parameter values
      for optionsdict in slicedoptionslist:
        # take a copy of the simulation options dictionary and set the 
        # values to the subset we've just extracted from the parameter sweep
        localsimoptionsdict = copy.deepcopy(simoptionsdict)
        localsimoptionsdict["values"] = optionsdict
        # append the simulation to the list
        self.simulations.append(simoptionsdict["type"](simulationpath, localsimoptionsdict, currentdirectory, tfdirectory, \
                                               logprefix=self.logprefix))

    self.runs = {}
    self.builds = {}
    self.collate()

    self.threadruns = []
    self.threadbuilds = []

    self.tests = [Test(name, code) for name, code in tests.iteritems()]

  def sliceoptionsdict(self, slicedoptionslist, diroptionsdict, optionsdict=None, i=0):
    '''Recursively slice the global options dictionary into a list of simulation optionsdictionaries 
       that contain a single parameter set from the full parameter sweep.'''

    # if this is the first level then set up an ordered dictionary 
    # (this gets repeatedly reused but we take deep copies into the list to ensure there's no overwriting)
    if optionsdict is None: optionsdict = collections.OrderedDict()

    # if we haven't reached the end yet...
    if i < len(diroptionsdict):
      # fetch the next sorted item
      if isinstance(diroptionsdict, collections.OrderedDict):
        item = diroptionsdict.items()[i]
      else:
        item = sorted(diroptionsdict.items())[i]
      # loop over the values provided in the list associated with item
      for strvalue in item[1]:
        # set the optionsdict to one of the values
        optionsdict[item[0]] = strvalue
        # and recurse to the next level of the global options dictionary
        self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1)
    # we've reached the lowest level of the diroptionsdict
    else:
      # create a deep copy of the optionsdict and append it to the list
      slicedoptionslist.append(copy.deepcopy(optionsdict))

  def collate(self):
    '''Collate the simulations based on their common run and build directories.'''

    # set up a dictionary of the runs
    self.runs = {}

    s = 0 # set up a counter
    # loop until we reach the end of the simulations list (which may at this point contain duplicates)
    while s < len(self.simulations):
      # if the rundirectory (which is used to identify unique simulations) is already registered then delete this simulation
      if self.simulations[s].rundirectory in self.runs:
        tmpsim = self.simulations.pop(s)
        del tmpsim
      # otherwise register the simulation in the run dictionary using its rundirectory as an identifier
      else:
        self.runs[self.simulations[s].rundirectory] = {'level':0, 'simulation':self.simulations[s]}
        # also register the simulation's dependencies at lower levels and get their deferment level
        dlevel = self.collatedependencies(self.simulations[s])
        # if this simulation is deferred then increment the dlevel
        if self.simulations[s].deferred: dlevel += 1
        # return the dlevel
        self.runs[self.simulations[s].rundirectory]['dlevel'] = dlevel
        # and increment the counter
        s += 1

    # set up a dictionary of the builds (empty if there are no builddirectories)
    self.builds = {}
    for run in self.runs.itervalues():
      if hasattr(run['simulation'], "builddirectory"):
        # this uses the builddirectory as the uniqueness identifier for simulations
        # all other information is copied from the runs dictionary
        self.builds[run['simulation'].builddirectory] = {'level':run['level'], 'dlevel':run['dlevel'],\
                                                         'simulation':run['simulation']}

  def collatedependencies(self, simulation, level=1):
    '''Collate the dependency simulations based on their common run directories.'''

    # set up the deferred level counter (0 is not deferred but any dependent of a deferred simulation must also be deferred)
    dlevel = 0
    # loop over the dependencies
    for s in xrange(len(simulation.dependencies)):
      # if the dependency rundirectory is already in the runs dictionary then we need to delete it and check its (d)level
      if simulation.dependencies[s].rundirectory in self.runs:
        # get the previous definition
        oldsim = self.runs[simulation.dependencies[s].rundirectory]
        # pop the new definition
        tmpsim = simulation.dependencies.pop(s)
        # add dependents of the new definition to the old definition
        oldsim['simulation'].dependents += tmpsim.dependents
        # delete the new (repeated definition)
        del tmpsim
        # insert the old definition into the slot previously occupied by the new one
        simulation.dependencies.insert(s, oldsim['simulation'])
        # ensure that the max level and dlevel are now used so that the dependency is run in an appropriate place
        oldsim['level'] = max(level, oldsim['level'])
        dlevel = max(oldsim['dlevel'], dlevel)
      # otherwise add the dependency to the runs dictionary
      else:
        self.runs[simulation.dependencies[s].rundirectory] = {'level':level, \
                                                              'simulation':simulation.dependencies[s]}
        # recursively call this function to collate the dependencies of this dependency
        ddlevel = self.collatedependencies(simulation.dependencies[s], level=level+1)
        # if the simulation is deferred then increment the deferred level counter from the dependencies
        if simulation.dependencies[s].deferred: ddlevel += 1
        # otherwise set the deferred level to that of the dependencies
        self.runs[simulation.dependencies[s].rundirectory]['dlevel'] = ddlevel
        # and update the deferred level we will be returning so we return the max of all the dependencies of this simulation
        dlevel = max(ddlevel, dlevel)
    # return the deferment level of the dependencies
    return dlevel

  def threadconfigure(self, queue, force=False):
    error = None
    for simulation in self.threadbuilds: 
      try:
        simulation.configure(force=force)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def configure(self, level=-1, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator(self.simulationselector(self.builds, level=level, dlevel=dlevel, types=types))
    for i in xrange(self.nthreads):
      queue = Queue.Queue()
      threadlist.append([threading.Thread(target=self.threadconfigure, args=[queue], kwargs={'force':force}), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)\n%s' % (ex_value.message, tb_str)
        raise ex_type(message) 

  def threadbuild(self, queue, force=False):
    error = None
    for simulation in self.threadbuilds: 
      try:
        simulation.build(force=force)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def build(self, level=-1, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator(self.simulationselector(self.builds, level=level, dlevel=dlevel, types=types))
    for i in xrange(self.nthreads):
      queue = Queue.Queue()
      threadlist.append([threading.Thread(target=self.threadbuild, args=[queue], kwargs={'force':force}), queue]) 
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)\n%s' % (ex_value.message, tb_str)
        raise ex_type(message) 

  def threadrun(self, queue, force=False):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.run(force=force)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def run(self, level=-1, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types))
    for i in xrange(self.nthreads):
      queue = Queue.Queue()
      threadlist.append([threading.Thread(target=self.threadrun, args=[queue], kwargs={'force':force}), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)\n%s' % (ex_value.message, tb_str)
        raise ex_type(message) 

    dlevel += 1
    if len(self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)) > 0:
      self.writeoptions(level=level, dlevel=dlevel, types=types)
      self.configure(level=level, dlevel=dlevel, types=types, force=force)
      self.build(level=level, dlevel=dlevel, types=types, force=force)
      self.run(level=level, dlevel=dlevel, types=types, force=force)

  def threadcheckpointrun(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.checkpointrun()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def checkpointrun(self, level=-1, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in xrange(self.nthreads):
      queue = Queue.Queue()
      threadlist.append([threading.Thread(target=self.threadcheckpointrun), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)\n%s' % (ex_value.message, tb_str)
        raise ex_type(message) 

  def listinput(self, level=-1, dlevel=0, types=None):
    runs = self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)[::-1]
    inputlist = set()
    for simulation in runs: inputlist.update(set(simulation.listinput()))
    return list(inputlist)

  def writeoptions(self, level=-1, dlevel=0, types=None):
    runs = self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)[::-1]
    for simulation in runs: simulation.writeoptions()  # not thread safe so don't thread it!

  def threadclean(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.clean()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def clean(self, level=-1, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in xrange(self.nthreads):
      queue = Queue.Queue()
      threadlist.append([threading.Thread(target=self.threadclean, args=[queue]), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)\n%s' % (ex_value.message, tb_str)
        raise ex_type(message) 

  def evaluatevariables(self, level=-1, dlevel=None, types=None):
    self.log("Assigning variables:")
    runs = self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)[::-1]
    varsdict = {}
    pathdict = {}
    varerror = False
    for simulation in runs: 
      try:
        simvarsdict = simulation.evaluatevariables()  # not thread safe so don't thread it!
      except SimulationsErrorVariable:
        varerror = True
        continue

      for name, value in simvarsdict.iteritems():
        # if we're not doing a parameter sweep then the values dictionary should be empty and we only accept
        # unique definitions of variable names
        if len(simulation.optionsdict["values"]) == 0:
          if name in varsdict:
            self.log("ERROR: Variable %s defined by simulation %s has been previously defined."%(name, simulation.name))
            varerror = True
            continue
          varsdict[name] = value
        # otherwise we are performing a parameter sweep and by definition there will be multiple values returned
        else:
          # check if we've come across this variable name before
          if name in varsdict:
            # the name being previously defined is only allowed if the simulation path is in the pathdict under this name
            prevdef = False
            try:
              prevdef = not pathdict[name]==simulation.path
            except KeyError:
              # if there is no entry in pathdict then this variable was previously defined by a single parameter simulation
              prevdef = True
            # throw an error if it has been previously defined by a different simulation (from a different parameter sweep)
            if prevdef:
              self.log("ERROR: Variable %s defined by simulation %s has been previously defined."%(name, simulation.name))
              varerror = True
              continue

          # otherwise if the name isn't in the varsdict then set up a nestedlist
          # and record the path in the pathdict so we can test for multiply defined variables later
          if not name in varsdict:
            varsdict[name] = NestedList(self.globaloptionsdict[simulation.path]["values"])
            pathdict[name] = simulation.path
          # set the values in the parameter sweep using the simulation options dictionary to index the nested list
          varsdict[name][simulation.optionsdict["values"]] = value

    if varerror: raise SimulationsErrorVariable

    return varsdict

  def runtests(self, force=False):
    '''Assign variables and run tests for all simulations.'''
    # set up a local variables dictionary
    varsdict = self.evaluatevariables()

    self.log("Running tests:")
    teststatus = []
    # change directory to the base directory of this group of simulations
    os.chdir(self.basedirectory)
    # loop over the tests
    for test in self.tests:
      self.log("Running %s:" % test.name)
      # run the test
      status = test.run(varsdict)
      if status == True:
        self.log("success.")
        teststatus.append('P')
      elif status == False:
        self.log("failure.")
        teststatus.append('F')
      else:
        self.log("failure (info == %s)." % status)
        teststatus.append('F')
    # change the directory back to the current directory
    os.chdir(self.currentdirectory)

    # print the test status report
    self.log(''.join(teststatus))

    # return the test status as a string
    return ''.join(teststatus)

  def test(self, force=False):
    '''Collates results from all tests in the SimulationBatch.'''
    
    # run the tests
    teststatus = self.runtests(force=force)

    # count the passes and failures
    passcount = teststatus.count('P')
    failcount = teststatus.count('F')
        
    # report
    if passcount + failcount > 0:
      self.log("Passes:   %d" % passcount)
      self.log("Failures: %d" % failcount)
    
    # exit with an error if any failures
    if failcount > 0:
      self.log("ERROR: at least one failure encountered while testing")
      raise SimulationsErrorTest

  def log(self, string):
    if self.logprefix is not None: string = self.logprefix+": "+string
    print string

  def simulationselector(self, simdict, level=-1, dlevel=None, types=None):
    return [value['simulation'] \
            for value in sorted(simdict.values(), key=lambda value: value['level']) \
            if value['level'] >= level and \
               (dlevel is None or value['dlevel'] == dlevel) and \
               (types is None or isinstance(value['simulation'], tuple(types)))]

#########################################################################################

class SimulationHarnessBatch(SimulationBatch):
  '''A derived SimulationBatch that takes in a list of spud based harnessfiles and sets up
     subgroups of SimulationBatches based on them.'''

  def __init__(self, harnessfiles, filename, currentdirectory, tfdirectory, nthreads=1):
    '''Initialize a SimulationHarnessBatch.'''

    # record the current directory from where the script calling this is being run
    self.currentdirectory = currentdirectory
    # the filename of the script calling this instance of this class
    self.filename = os.path.basename(filename)
    # work out the directory of the file (used for running the tests in)
    dirname = os.path.dirname(filename)
    if os.path.isabs(dirname):
      self.basedirectory = os.path.normpath(dirname)
    else:
      self.basedirectory = os.path.normpath(os.path.join(currentdirectory, dirname))
    # no log prefix at this level
    self.logprefix = None
    # number of threads we are to run processes over
    self.nthreads = nthreads

    # set up a list sub SimulationBatches - one for each harnessfile so that tests can be grouped together
    self.simulationtestgroups = []
    # establish the global options dictionary for this batch
    self.globaloptionsdict = {}
    # loop over the list of harnessfiles that has been passed in.
    for harnessfile in harnessfiles:
      # set up an optionsdictionary for this harnessfile
      harnessfileoptionsdict = {}
      dirname = os.path.dirname(harnessfile)
      # the harnessfile should be a .shml libspud compatible file
      libspud.load_options(harnessfile)
      # loop over the simulations
      for s in xrange(libspud.option_count("/simulations/simulation")):
        simulation_optionpath = "/simulations/simulation["+`s`+"]"
        simulation_name = libspud.get_option(simulation_optionpath+"/name")
        
        # name of options file and normalize path to it relative to the path to the harnessfile itself
        simulation_input_file = libspud.get_option(simulation_optionpath+"/input_file")
        if os.path.isabs(simulation_input_file):
          simulation_path = os.path.normpath(simulation_input_file)
        else:
          simulation_path = os.path.normpath(os.path.join(dirname, simulation_input_file))

        # fetch required input and output for this simulation (used to copy files into rundirectory 
        # and establish if run has successfully completed
        required_input  = self.getrequiredfiles(simulation_optionpath+"/required_input", dirname)
        required_output = self.getrequiredfiles(simulation_optionpath+"/required_output")

        # get the dictionary of the parameters we are sweeping over in this simulation and the number of times each run should
        # be performed (e.g. for benchmarking purposes)
        parameter_values, parameter_updates, parameter_builds = self.getparameters(simulation_optionpath+"/parameter_sweep")
        try:
          nruns = libspud.get_option(simulation_optionpath+"/parameter_sweep/number_runs")
        except libspud.SpudKeyError:
          nruns = 1

        # get the parameters for any checkpoint pickups
        checkpoint_values, checkpoint_updates, checkpoint_builds = self.getparameters(simulation_optionpath+"/checkpoint")
        # don't do anything with checkpoint_builds

        # get the details of variables we want to examine from this simulation
        variables = self.getvariables(simulation_optionpath+"/variables")

        # add details of this simulation to the harness file options dictionary
        harnessfileoptionsdict[simulation_path] = {}
        harnessfileoptionsdict[simulation_path]["name"]    = simulation_name
        harnessfileoptionsdict[simulation_path]["type"]    = Simulation
        harnessfileoptionsdict[simulation_path]["input"]   = required_input
        harnessfileoptionsdict[simulation_path]["output"]  = required_output
        harnessfileoptionsdict[simulation_path]["values"]  = parameter_values
        harnessfileoptionsdict[simulation_path]["updates"] = parameter_updates
        harnessfileoptionsdict[simulation_path]["builds"]  = parameter_builds
        harnessfileoptionsdict[simulation_path]["nruns"]   = nruns
        harnessfileoptionsdict[simulation_path]["checkpoint_values"]  = checkpoint_values
        harnessfileoptionsdict[simulation_path]["checkpoint_updates"] = checkpoint_updates
        harnessfileoptionsdict[simulation_path]["variables"] = variables
        # if the simulation has dependencies then recursively get their options too
        if libspud.have_option(simulation_optionpath+"/dependencies"):
          harnessfileoptionsdict[simulation_path]["dependencies"] = \
             self.getdependencies(simulation_optionpath+"/dependencies", dirname, parameter_values.keys())

      # fetch the tests that are included in this harness file
      tests = self.gettests("/tests")

      # append the details of this harnessfile to the list of simulation test groups
      self.simulationtestgroups.append(SimulationBatch(harnessfileoptionsdict, \
                                                       harnessfile, currentdirectory, tfdirectory, \
                                                       tests=tests, nthreads=nthreads))
      # add to the global options dictionary for this 'super' SimulationBatch
      self.globaloptionsdict.update(harnessfileoptionsdict)
      # clear the options so the next shml file can be loaded
      libspud.clear_options()

    # set up a list of simulations from all the sub SimulationBatches we just set up
    self.simulations = []
    for simulationtestgroup in self.simulationtestgroups:
      for simulation in simulationtestgroup.simulations:
        self.simulations.append(simulation)

    # from that list of simulations collate all the unique ones in dictionaries of runs and builds
    self.runs = {}
    self.builds = {}
    self.collate()

    # some final lists that are necessary for running, building and testing
    self.threadruns = []
    self.threadbuilds = []
    self.tests = []

  def test(self, force=False):
    '''Collates results from all tests in the sub SimulationBatches.'''

    teststatus = ''
    failinggroups = []
    varerror = False
    # loop over the simulation groups
    for group in self.simulationtestgroups:
      # run the tests on each group
      try:
        groupstatus = group.runtests(force=force)
      except SimulationsErrorVariable:
        varerror = True
        continue
      # append any failing groups to the summary list
      if groupstatus.count('F') > 0:
        failinggroups.append([group.logprefix, groupstatus])
      # append the status to the collated list
      teststatus = teststatus + groupstatus

    passcount = teststatus.count('P') # count the passes
    failcount = teststatus.count('F') # count the failures
    
    # if there were failures then summarize the failures
    if failcount > 0:
       print
       print "Summary of test problems with failures:"
       for group in failinggroups:
          print group[0]+':', group[1]
       print
    
    # report the overall number of passes and failures
    if passcount + failcount > 0:
      self.log("Passes:   %d" % passcount)
      self.log("Failures: %d" % failcount)
    
    # exit with an error if any failed
    if varerror: raise SimulationsErrorVariable
    if failcount > 0: raise SimulationsErrorTest

  def getrequiredfiles(self, optionpath, dirname=None):
     '''Get a list of the required files listed under the optionpath.  Run python scripts in dirname if supplied.'''

     # regex to split up filename lists
     # accept comma, semicolon or space delimiters
     r = re.compile(r'(?:[^,; ])+')

     required_files = []
     # loop over the filenames
     for f in xrange(libspud.option_count(optionpath+"/filenames")):
       filename_optionpath = optionpath+"/filenames["+`f`+"]"
       # try getting the files as a string
       try:
         filenames = libspud.get_option(filename_optionpath+"/string")
         required_files += r.findall(filenames)
       # otherwise run a python string to get the filename(s)
       except libspud.SpudKeyError:
         filename_name = libspud.get_option(filename_optionpath+"/name")
         filename_code = libspud.get_option(filename_optionpath+"/python")
         var = Variable(filename_name, filename_code)
         varsdict = {}
         if dirname is not None: os.chdir(dirname)
         var.run(varsdict)
         os.chdir(self.currentdirectory)
         filenames = varsdict[filename_name]
         # we accept either a string or a list of filenames from python
         if isinstance(filenames, str):
           required_files.append(filenames)
         else:
           required_files += filenames

     # if dirname has been supplied then make the filenames relative to it
     if dirname is not None:
       for f in xrange(len(required_files)):
         filename = required_files[f]
         if os.path.isabs(filename):
           required_files[f] = os.path.normpath(filename)
         else:
           required_files[f] = os.path.normpath(os.path.join(dirname, filename))

     # return the required files
     return required_files

  def getcommands(self, optionpath):
     commands = []
     for c in xrange(libspud.option_count(optionpath+"/command")):
       command_optionpath = optionpath+"/command["+`c`+"]"
       commands.append(libspud.get_option(command_optionpath).split())
     return commands

  def getvariables(self, optionpath):
     variables  = collections.OrderedDict()
     for v in xrange(libspud.option_count(optionpath+"/variable")):
       variable_optionpath = optionpath+"/variable["+`v`+"]"
       variable_name = libspud.get_option(variable_optionpath+"/name")
       variables[variable_name] = libspud.get_option(variable_optionpath)
     return variables

  def gettests(self, optionpath):
     tests  = collections.OrderedDict()
     for t in xrange(libspud.option_count(optionpath+"/test")):
       test_optionpath = optionpath+"/test["+`t`+"]"
       test_name = libspud.get_option(test_optionpath+"/name")
       tests[test_name] = libspud.get_option(test_optionpath)
     return tests

  def getparameters(self, optionpath, parent_parameters=None):
     '''Fetch parameters from /parameters child of the given optionpath.  
        Test that the parameters are available in the parent if parent_parameters present.'''

     # regular expression to split up values string into a list
     # the list may be comma(,), semicolon(;), or space ( ) delimited
     # if the list is of bracketed items (e.g. tuples) then any delimiters within 
     # the brackets are preserved
     # brackets may be (), [] or {}
     #NODE                     EXPLANATION
     #--------------------------------------------------------------------------------
     #  (?:                      group, but do not capture (1 or more times
     #                           (matching the most amount possible)):
     #--------------------------------------------------------------------------------
     #    [^,; ([{]                any character except: ',', ';', ' ',
     #                             '(', '[', '{'
     #--------------------------------------------------------------------------------
     #   |                        OR
     #--------------------------------------------------------------------------------
     #    \(                       '('
     #--------------------------------------------------------------------------------
     #    [^)]*                    any character except: ')' (0 or more
     #                             times (matching the most amount
     #                             possible))
     #--------------------------------------------------------------------------------
     #    \)                       ')'
     #--------------------------------------------------------------------------------
     #   |                        OR
     #--------------------------------------------------------------------------------
     #    \[                       '['
     #--------------------------------------------------------------------------------
     #    [^]]*                    any character except: ']' (0 or more
     #                             times (matching the most amount
     #                             possible))
     #--------------------------------------------------------------------------------
     #    \]                       ']'
     #--------------------------------------------------------------------------------
     #   |                        OR
     #--------------------------------------------------------------------------------
     #    \{                       '{'
     #--------------------------------------------------------------------------------
     #    [^}]*                    any character except: '}' (0 or more
     #                             times (matching the most amount
     #                             possible))
     #--------------------------------------------------------------------------------
     #    \}                       '}'
     #--------------------------------------------------------------------------------
     #  )+                       end of grouping
     # - from http://rick.measham.id.au/paste/explain.pl
     r = re.compile(r'(?:[^,; ([{]|\([^)]*\)|\[[^]]*\]|\{[^}]*\})+')
    
     # set up ordered dictionaries to remember the order in which parameters are specified
     parameter_values  = collections.OrderedDict()
     parameter_updates = collections.OrderedDict()
     parameter_builds  = collections.OrderedDict()
     # loop over all the parameters
     for p in xrange(libspud.option_count(optionpath+"/parameter")):
       parameter_optionpath = optionpath+"/parameter["+`p`+"]"
       parameter_name = libspud.get_option(parameter_optionpath+"/name")
       
       # throw an error if the parameter is not in the parent list (if present)
       if parent_parameters is not None:
         if parameter_name not in parent_parameters:
           self.log("ERROR: dependency parameter not in parent simulation:")
           self.log("parameter %(parameter_name)s (%(parameter_optionpath)s)"%\
                 {"parameter_name": parameter_name, "parameter_optionpath": parameter_optionpath})
           
           raise SimulationsErrorInitialization
       
       # try to get a list of the parameter values (not necessarily specified for child parameters so except failure)
       try:
         parameter_values[parameter_name]  = r.findall(libspud.get_option(parameter_optionpath+"/values"))
       except libspud.SpudKeyError:
         pass # insert nothing if we haven't specified any values

       # try to get a list of the parameter update instructions (not necessarily specified for parent parameters so except failure)
       try:
         parameter_updates[parameter_name] = libspud.get_option(parameter_optionpath+"/update")
         parameter_builds[parameter_name]  = libspud.have_option(parameter_optionpath+"/update/single_build")
       except libspud.SpudKeyError:
         parameter_updates[parameter_name] = None # insert None if we haven't said how to update
         parameter_builds[parameter_name]  = True # not updating so we don't need multiple builds

     # return the OrderedDict objects for the values and updates
     return parameter_values, parameter_updates, parameter_builds

  def getdependencies(self, optionpath, dirname, parent_parameters):
     dependencies_options = {}

     for d in xrange(libspud.option_count(optionpath+"/simulation")):
        simulation_optionpath = optionpath+"/simulation["+`d`+"]"
        dependencies_options.update(self.getdependency(simulation_optionpath, dirname, parent_parameters))
     
     for d in xrange(libspud.option_count(optionpath+"/run")):
        run_optionpath = optionpath+"/run["+`d`+"]"
        dependencies_options.update(self.getdependency(run_optionpath, dirname, parent_parameters, run=True))
     
     return dependencies_options

  def getdependency(self, optionpath, dirname, parent_parameters, run=False):
     dependency_options = {}
     name = libspud.get_option(optionpath+"/name")

     input_file = libspud.get_option(optionpath+"/input_file")
     if os.path.isabs(input_file):
       path = os.path.normpath(input_file)
     else:
       path = os.path.normpath(os.path.join(dirname, input_file))
     
     required_input  = self.getrequiredfiles(optionpath+"/required_input", dirname)
     required_output = self.getrequiredfiles(optionpath+"/required_output")

     parameter_values, parameter_updates, parameter_builds = self.getparameters(optionpath+"/parameter_sweep", \
                                                                                parent_parameters=parent_parameters)
     # we do nothing with parameter_values... it should be empty
     try:
       nruns = libspud.get_option(optionpath+"/parameter_sweep/number_runs")
     except libspud.SpudKeyError:
       nruns = 1
     
     variables = self.getvariables(optionpath+"/variables")

     dependency_options[path] = {}
     dependency_options[path]["name"]      = name
     if run:
       dependency_options[path]["type"]    = Run
     else:
       dependency_options[path]["type"]    = Simulation
     dependency_options[path]["input"]     = required_input
     dependency_options[path]["output"]    = required_output
     dependency_options[path]["updates"]   = parameter_updates
     dependency_options[path]["builds"]    = parameter_builds
     dependency_options[path]["nruns"] = nruns
     dependency_options[path]["variables"] = variables
     if libspud.have_option(optionpath+"/dependencies"):
       dependency_options[path]["dependencies"] = \
                           self.getdependencies(optionpath+"/dependencies", dirname, \
                                                parameter_updates.keys())

     if run:
       dependency_options[path]["run"] = self.getcommands(optionpath+"/commands")

     return dependency_options
