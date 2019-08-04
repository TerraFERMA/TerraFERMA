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
import gzip
import shutil
import threading
from queue import Queue
import copy
import glob
import re
import collections
from buckettools.threadlibspud import *
import traceback
from functools import reduce
import operator
from string import Template as template

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

class TestOrVariableException(Exception):
  def __init__(self, msg = "Error while evaluating test or variable."):
    self.message = msg

####################################################################################

class ThreadIterator(list):
  '''A thread-safe iterator over a list.'''
  def __init__(self, seq):
    self.list=list(seq)
    self.lock=threading.Lock()

  def __iter__(self):
    return self

  def __next__(self):
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
        tmpdict = copy.deepcopy(varsdict) # don't let the test code modify the variables
        try:
          exec(self.code, tmpdict)
        except AssertionError:
          # in case of an AssertionError, we assume the test has just failed
          return False
        except:
          message = "Test computation raised an exception."+os.linesep
          message += "-" * 80 + os.linesep
          for (lineno, line) in enumerate(self.code.splitlines()):
            message + "%4d  %s" % (lineno+1, line) + os.linesep
          message += "-" * 80 + os.linesep
          message += traceback.format_exc()
          message += "-" * 80
          raise TestOrVariableException(msg=message)
        else:
          return True

####################################################################################

class Variable(TestOrVariable):
    """A variable definition for use in tests"""
    def run(self, varsdict):
        try:
          exec(self.code, varsdict)
        except:
          message = "Variable computation raised an exception."+os.linesep
          message += "-" * 80 + os.linesep
          for (lineno, line) in enumerate(self.code.splitlines()):
            message += "%4d  %s" % (lineno+1, line) + os.linesep
          message += "-" * 80 + os.linesep
          message += traceback.format_exc()
          message += "-" * 80 + os.linesep
          raise TestOrVariableException(msg=message)
        else:
          if self.name not in list(varsdict.keys()):
            message = "Variable (%s) not defined by computation."%(self.name)+os.linesep
            message += "-" * 80 + os.linesep
            for (lineno, line) in enumerate(self.code.splitlines()):
              message += "%4d  %s" % (lineno+1, line) + os.linesep
            message += "-" * 80 + os.linesep
            raise TestOrVariableException(msg=message)

####################################################################################

class NestedList(list):
  '''A class that implements a nested list structure based on a dictionary of parameters.  
     The resulting nested list can be indexed by a dictionary of parameters that works out which
     indices to index into the nested list structure.'''

  def __init__(self, parameters=None):
    '''Initialize a nested list.'''
    # initialize the base list class
    super(NestedList, self).__init__()

    # record the parameter dictionary or set up an empty one
    if parameters is None: parameters = collections.OrderedDict()
    self.parameters = parameters

    # if the parameters have been provided in an ordered dictionary then
    # we respect that order, otherwise we sort the keys alphabetically:
    if isinstance(self.parameters, collections.OrderedDict):
      self.sortedkeys = list(self.parameters.keys())
    else:
      self.sortedkeys = sorted(self.parameters.keys())

    # create the nested list structure, beneath the top level list object, self:
    self.createnestedlist(self)

  def createnestedlist(self, parentlist, level=0):
    '''Recursively set up a nested list.'''
    # if we've reached the bottom of the dictionary then append Nones and done recurse
    if level == len(self.sortedkeys)-1:
      for val in self.parameters[self.sortedkeys[level]]: parentlist.append(None)
    # if we haven't reached the bottom then loop over the values and recursively call
    # this function to set up the sublevels of the nested list
    elif level < len(self.sortedkeys)-1:
      for val in self.parameters[self.sortedkeys[level]]: 
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
      paramvals = self.parameters[key]
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
        indicies = list(range(len(paramvals)))
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
  def __init__(self, path, optionsdict, currentdirectory, tfdirectory, batchdirectory, \
               logprefix=None, dependents=[]):
    self.optionsdict = optionsdict

    filename = os.path.basename(path)
    self.filename, self.ext = os.path.splitext(filename)
    self.inputfilename = self.filename
    self.inputext = self.ext
    self.spudfile = False
    if "spudfile" in self.optionsdict:
      self.spudfile = self.optionsdict["spudfile"]

    self.name = self.optionsdict["name"]

    self.logprefix = logprefix

    self.currentdirectory = currentdirectory
    self.batchdirectory   = batchdirectory

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

    self.variables = [Variable(name, code) for name, code in self.optionsdict["variables"].items()]

    self.alreadyrun = False

    self.dependencies = []
    if "dependencies" in self.optionsdict:
       for dependencypath, depoptionsdict in self.optionsdict["dependencies"].items():
         dependencyoptions = self.getdependencyoptions(depoptionsdict["procscales"])
         depoptionsdict["values"] = dependencyoptions["values"]
         depoptionsdict["procscales"] = dependencyoptions["procscales"]
         depoptionsdict["value_indices"] = dependencyoptions["value_indices"]
         depoptionsdict["value_lengths"] = dependencyoptions["value_lengths"]
         self.dependencies.append(depoptionsdict["type"](dependencypath, depoptionsdict, \
                                                         currentdirectory, tfdirectory, batchdirectory, \
                                                         logprefix=logprefix, dependents=[self]))

    # if the input options file is an the output of a dependency then we have to
    # defer generating, configuring and building until runtime
    filepaths = [depoutput_k  \
                 for depoutput_k, depoutput_v in list(self.getdependencyrequiredoutput().items()) \
                 if self.filename+self.ext == os.path.basename(depoutput_v)]
    self.deferred = len(filepaths)>0
    if self.deferred:
      filepath = filepaths[0]
      dirname = os.path.dirname(filepath)
      self.inputfilename, self.inputext = os.path.splitext(os.path.basename(filepath))
      if os.path.isabs(dirname):
        self.runinputdirectory = os.path.normpath(dirname)
      else:
        self.runinputdirectory = os.path.normpath(os.path.join(currentdirectory, dirname))

  def writeoptions(self):

    try:
      os.makedirs(self.rundirectory)
    except OSError:
      pass
      
    if self.spudfile:
      threadlibspud.load_options(os.path.join(self.runinputdirectory, self.inputfilename+self.inputext))

      self.updateoptions()

      libspud.write_options(os.path.join(self.rundirectory, self.filename+self.ext))
      threadlibspud.clear_options()
    else:
      inputfilepath = os.path.join(self.runinputdirectory, self.inputfilename+self.inputext)
      gz = False
      try:
        inputfile = open(inputfilepath)
        inputstr = inputfile.read()
      except UnicodeDecodeError:
        inputfile = gzip.open(inputfilepath, 'rt', encoding='utf-8')
        inputstr = inputfile.read()
        gz = True
      inputfile.close()
      
      valuesdict = {"input_file":inputstr}
      self.updateoptions(valuesdict=valuesdict)

      inputstr = valuesdict["input_file"]
      outputfilename = os.path.join(self.rundirectory, self.filename+self.ext)
      if gz:
        outputfile = gzip.open(outputfilename, 'wt')
      else:
        outputfile = open(outputfilename, 'w')
      outputfile.write(inputstr)
      outputfile.close()

  def updateoptions(self, valuesdict=None, prefix=""):
    os.chdir(self.basedirectory)
    sys.path.append(self.basedirectory)
    valuesdict = self.getvaluesdict(valuesdict=valuesdict, prefix=prefix)
    for paramname, update in self.optionsdict[prefix+"updates"].items():
      if update is not None:
        try:
          exec(update, valuesdict)
        except libspud.SpudNewKeyWarning as e:
          self.log("ERROR: spud raised a new key warning:")
          self.log("%s"%e)
          self.log("on parameter: %s"%paramname)
          self.log("with update:"+os.linesep+"%s"%update) 
          raise SimulationsErrorWriteOptions
        except Exception as e:
          self.log("ERROR: evaluating parameters raised an exception:")
          self.log("%s"%e)
          self.log("on parameter: %s"%paramname)
          self.log("with update:"+os.linesep+"%s"%update) 
          raise SimulationsErrorWriteOptions
    sys.path.remove(self.basedirectory)
    os.chdir(self.currentdirectory)

  def getvaluesdict(self, valuesdict=None, prefix=""):
    if valuesdict is None: valuesdict={}
    valuesdict["input_filename"] = os.path.join(self.runinputdirectory, self.inputfilename+self.inputext)
    valuesdict["_self"] = self
    for k,v in self.optionsdict[prefix+"values"].items():
      if k in valuesdict:
        self.log("ERROR: in getvaluesdict, %s multiply defined"%(k))
        if k=="input_file":
          self.log("parameter cannot be named input_file")
        if k=="input_filename":
          self.log("parameter cannot be named input_filename")
        if k=="_self":
          self.log("parameter cannot be named _self")
        raise SimulationsErrorWriteOptions
      valuesdict[k] = v
    return valuesdict

  def resolvepythonrequiredfiles(self, required_files_python, dirname=None):
     '''Get a list of files from the supplied list of filenames and python code.'''

     required_files = {}
     for filename_name, filename_code in required_files_python:
       var = Variable(filename_name, filename_code)
       varsdict = self.getvaluesdict()
       os.chdir(self.batchdirectory)
       try:
         var.run(varsdict)
       except TestOrVariableException as exc:
         self.log("ERROR: error while evaluating filename(s) %s."%(filename_name))
         self.log(exc.message)
         sys.exit(1)
       os.chdir(self.currentdirectory)
       filenames = varsdict[filename_name]
       # we accept either a string or a list of filenames from python
       if isinstance(filenames, str):
         required_files[filenames] = filenames
       elif isinstance(filenames, list):
         required_files.update({filename : filename for filename in filenames})
       elif isinstance(filenames, dict):
         required_files.update(filenames)
       else:
         self.log("ERROR: Unknown format of required input supplied.")
         raise SimulationsErrorInitialization
         

     # if dirname has been supplied then make the filenames relative to it
     if dirname is not None:
       new_required_files = {(os.path.normpath(filename_k) if os.path.isabs(filename_k) 
                              else os.path.normpath(os.path.join(dirname, filename_k))):filename_v
                              for filename_k, filename_v in list(required_files.items())}
       required_files = new_required_files

     # return the required files
     return required_files

  def configure(self, force=False):
    return None

  def build(self, force=False):
    return None

  def run(self, force=False, rundependencies=True):
    
    self.lock.acquire()

    error = False
    
    if rundependencies: 
      for dependency in self.dependencies: dependency.run(force=force)

    self.log("Checking in directory: %s"%(os.path.relpath(self.rundirectory, self.currentdirectory)))
    if not self.alreadyrun and (not self.optionsdict["run_when"]["never"] or force):
      commands = self.getcommands()

      valuesdict=self.getvaluesdict()

      for r in range(self.nruns):
        requiredinput = self.getrequiredinput(r)
        requiredoutput = self.getrequiredoutput(r)

        dirname = os.path.join(self.rundirectory, "run_"+repr(r).zfill(len(repr(self.nruns))))
        try:
          os.makedirs(dirname)
        except OSError:
          pass

        input_changed = False
        if self.optionsdict["run_when"]["input_changed"]:
          for filepath_k, filepath_v in requiredinput.items():
            try:
              try:
                fh = open(os.path.join(dirname, os.path.basename(filepath_v))).read()
              except UnicodeDecodeError:
                fh = gzip.open(os.path.join(dirname, os.path.basename(filepath_v)), 'rt', encoding='utf-8').read()
              checksum = hashlib.md5(fh.encode('utf-8')).hexdigest()
            except:
              checksum = None
            try:
              try:
                fh = open(filepath_k).read()
              except UnicodeDecodeError:
                fh = gzip.open(filepath_k, 'rt', encoding='utf-8').read()
              input_changed = input_changed or checksum != hashlib.md5(fh.encode('utf-8')).hexdigest()
            except IOError:
              self.log("WARNING: Unable to open %s"%(filepath_k))
              input_changed = True

        output_missing = False
        if self.optionsdict["run_when"]["output_missing"]:
          for filepath_k, filename_v in requiredoutput.items():
            try:
              output_file = open(os.path.join(dirname, filepath_k))
              output_file.close()
            except IOError:
              output_missing = True
          if len(requiredoutput)==0: output_missing=True # don't know what output is needed so we have to force
          
        if output_missing or input_changed or force or self.optionsdict["run_when"]["always"]:
          self.log("  Running in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))
          # file has changed or a recompilation is necessary
          for filepath_k, filepath_v in requiredoutput.items():
            try:
              os.remove(os.path.join(dirname, filepath_k))
            except OSError:
              pass
          for filepath_k, filepath_v in requiredinput.items():
            try:
              shutil.copy(filepath_k, os.path.join(dirname, os.path.basename(filepath_v)))
            except IOError:
              self.log("WARNING: required input (%s) not found, continuing anyway."%(filepath_k))
          
          env = copy.deepcopy(os.environ)
          try:
            env["PYTHONPATH"] = ":".join([dirname, env["PYTHONPATH"]])
          except KeyError:
            env["PYTHONPATH"] = dirname
          env["PWD"] = dirname

          for i in range(len(commands)):
            command = commands[i]
            tcommand = [template(c).safe_substitute(valuesdict) for c in command]
            logf = '_'+self.filename+self.ext+'.'+repr(i)+'.log'
            retvalue = self.runcommand(tcommand, dirname, logfilename=logf)
            error = retvalue != 0
            if error: break
          
          if error:
            # There's been an error, append failed output
            for filepath_k, filename_v in requiredoutput.items():
              if os.path.isfile(os.path.join(dirname, filepath_k)):
                shutil.move(os.path.join(dirname, filepath_k), os.path.join(dirname, filepath_k+".fail"))
          else:
            # Cleanup previous errors (if any)
            for filepath_k, filename_v in requiredoutput.items():
              if os.path.isfile(os.path.join(dirname, filepath_k+".fail")):
                os.remove(os.path.join(dirname, filepath_k+".fail"))

          self.log("  Finished in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))

      self.alreadyrun = True

    self.lock.release()

    if error: raise SimulationsErrorRun

  def checkpointrun(self, index=-1):
    pass

  def cleanbuild(self):
    pass

  def cleanrun(self):
    try:
      shutil.rmtree(os.path.join(self.basedirectory, self.filename+self.ext+".run"))
    except OSError:
      pass

  def clean(self):
    self.cleanrun()

  def checkpointclean(self):
    for r in range(self.nruns):
      dirname = os.path.join(self.rundirectory, "run_"+repr(r).zfill(len(repr(self.nruns))), "checkpoint")
      try:
        shutil.rmtree(dirname)
      except OSError:
        pass

  def evaluatevariables(self):
    '''Evaluate the variables for this simulation and its dependencies.'''

    error = False
    # set up the variables dictionary
    varsdict = {}
    # loop over all the runs
    for r in range(self.nruns):
      # change to the run directory to run the python variable assignment
      dirname = os.path.join(self.rundirectory, "run_"+repr(r).zfill(len(repr(self.nruns))))
      try:
        os.chdir(dirname)
      except OSError:
        self.log("WARNING: could not change to directory %s"%(dirname))

      # loop over the variables
      for var in self.variables:
        if self.nruns > 1: varsdict = []
        # set up a temporary variable dictionary for this variable calculation
        tmpdict  = self.getvaluesdict()
        # try running the variable assignment
        try:
          var.run(tmpdict)
        except TestOrVariableException as exc:
          self.log("ERROR: failure while calculating variable %s." % (str(var.name)))
          self.log(exc.message)
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
    for line in string.splitlines():
      if self.logprefix is not None: line = self.logprefix+": "+line
      sys.stdout.write(line+os.linesep)
    sys.stdout.flush()

  def getrundirectory(self):
    '''Return the run directory for this simulation.'''

    # start at the base directory
    rundirectory = os.path.join(self.basedirectory, self.filename+self.ext+".run")

    # preserve the order of an ordered dictionary or use alphabetical
    if isinstance(self.optionsdict["values"], collections.OrderedDict):
      sorteditems = list(self.optionsdict["values"].items())
    else:
      sorteditems = sorted(self.optionsdict["values"].items())

    # append directory levels for all parameter values
    for directory in [item[0]+"_"+item[1] for item in sorteditems]: 
      rundirectory = os.path.join(rundirectory, directory)

    # return the run directory
    return rundirectory

  def getdependencyoptions(self, procscales):
    suboptionsdict = {"values" : collections.OrderedDict(), "procscales" : collections.OrderedDict(), 
                      "value_indices" : collections.OrderedDict(), "value_lengths" : collections.OrderedDict() }
    for option, procscalelist in procscales.items(): 
      suboptionsdict["values"][option] = self.optionsdict["values"][option]
      suboptionsdict["value_indices"][option] = self.optionsdict["value_indices"][option]
      suboptionsdict["value_lengths"][option] = self.optionsdict["value_lengths"][option]

      if procscalelist is not None:
        if len(procscalelist) != self.optionsdict["value_lengths"][option]:
          self.log("ERROR: number of process scales of dependency does not equal the number of values in root:")
          self.log("parameter %(parameter_name)s"%{"parameter_name": option})
          raise SimulationsErrorInitialization

        suboptionsdict["procscales"][option] = procscalelist[self.optionsdict["value_indices"][option]]
      else:
        suboptionsdict["procscales"][option] = 1

    return suboptionsdict

  def getrequiredoutput(self, run):
    requiredoutput = {}
    if "output" in self.optionsdict:
      outputfiles = self.optionsdict["output"]
      for outputfile in outputfiles: requiredoutput[outputfile] = outputfile
    if "output_python" in self.optionsdict:
      requiredoutput.update(self.resolvepythonrequiredfiles(self.optionsdict["output_python"]))
    return requiredoutput
 
  def getdependencyrequiredoutput(self, run=0):
    requiredoutput = {}
    for dependency in self.dependencies:
      if self.nruns == dependency.nruns: 
        rundir = "run_"+repr(run).zfill(len(repr(self.nruns)))
      else:
        rundir = "run_"+repr(0).zfill(len(repr(dependency.nruns)))
      dependencyoutput = dependency.getrequiredoutput(run)
      for filename_k, filename_v in list(dependencyoutput.items()):
        requiredoutput[os.path.join(dependency.rundirectory, rundir, filename_k)] = filename_v
      requiredoutput.update(dependency.getdependencyrequiredoutput(run))
    return requiredoutput

  def listinput(self):
    requiredinput = []
    if not self.deferred: requiredinput.append(os.path.join(self.basedirectory, self.filename+self.ext))
    if "input" in self.optionsdict:
      requiredinput += self.optionsdict["input"]
    if "input_python" in self.optionsdict:
      requiredinput += \
         list(self.resolvepythonrequiredfiles(self.optionsdict["input_python"], dirname=self.batchdirectory).keys())
    return requiredinput

  def getrequiredinput(self, run, build=False):
    input_filename = os.path.join(self.rundirectory, self.filename+self.ext)
    requiredinput = {input_filename:input_filename}
    inputkey = "input"
    if build: inputkey = "build_input"
    # get any input specified in the optionsdict
    if inputkey in self.optionsdict: 
      for filepath in self.optionsdict[inputkey]:
        # filter out any filenames that have already had their destination filename added (i.e. make sure our latest input file is used)
        if os.path.basename(filepath) not in [os.path.basename(inputpath) for inputpath in list(requiredinput.values())]:
          requiredinput[filepath] = filepath
    if inputkey+"_python" in self.optionsdict:
      for filepath_k, filepath_v in self.resolvepythonrequiredfiles(self.optionsdict[inputkey+"_python"], dirname=self.batchdirectory).items():
        # filter out any filenames that have already been added (i.e. make sure our latest input file is used)
        if os.path.basename(filepath_v) not in [os.path.basename(inputpath) for inputpath in list(requiredinput.values())]:
          requiredinput[filepath_k] = filepath_v
    if not build:
      # get any output from dependencies
      for filepath_k, filepath_v in self.getdependencyrequiredoutput(run).items():
        # filter out any filenames that have already been added (highest dependencies take priority over lower ones)
        if os.path.basename(filepath_v) not in [os.path.basename(inputpath) for inputpath in list(requiredinput.values())]:
          requiredinput[filepath_k] = filepath_v
    return requiredinput

  def getcommands(self):
    commands = []
    if "run" in self.optionsdict: commands = self.optionsdict["run"]
    return commands

  def runcommand(self, command, dirname, verbose=False, exception=None, logfilename=None, env=None):
    """Run a command through subprocess"""
    # start the command
    try:
      p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=dirname, env=env)
    except Exception as exc:
      self.log("  ERROR: %s raised an exception before running:"%(os.path.split(command[0])[-1]))
      self.log("         %s"%str(exc))
      self.log("         in directory: %s"%(dirname))
      self.log("         using command: %s"%(" ".join(command)))
      if exception is not None:
        raise exception
      else:
        return 1
    output = bytearray()
    # open a log
    if logfilename is not None:
      logfile = open(os.path.join(dirname, logfilename), "w")
    else:
      logfile = None
    # continually iterate over the lines of stdout and pipe them to the logs
    for line in iter(p.stdout.readline, b''):
      if verbose: self.log(" "*2+line.decode())
      if logfile is not None: logfile.write(line.decode())
      output.extend(line)
    if logfile is not None: logfile.close() # close the logfile
    retcode = p.wait()
    # if we've got a non-zero return code report the error and raise an exception or return a failure code
    if retcode:
      self.log("  ERROR: %s returned: %d"%(os.path.split(command[0])[-1], retcode))
      self.log("         in directory: %s"%(dirname))
      self.log("         using command: %s"%(" ".join(command)))
      for line in output.decode().splitlines():
        self.log(" "*2 + line + os.linesep)
      if exception is not None: 
        raise exception
      else:
        return 1
    else:
      # otherwise return a success signal
      return 0

################################################################################################

class Simulation(Run):
  def __init__(self, path, optionsdict, currentdirectory, tfdirectory, batchdirectory, \
               logprefix=None, dependents=[]):

    Run.__init__(self, path, optionsdict, currentdirectory, tfdirectory, batchdirectory, \
                 logprefix=logprefix, dependents=dependents)

    self.tfdirectory = tfdirectory
    self.builddirectory = self.getbuilddirectory()

  def writeoptions(self):
    
    threadlibspud.load_options(os.path.join(self.runinputdirectory, self.inputfilename+self.inputext))
    
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
      os.makedirs(os.path.join(basedir, "checkpoint"))
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

    # run twice as cmake seems to modify flags on second run
    self.log("Configuring in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))
    command = ["cmake", \
                        "-DOPTIONSFILE="+os.path.join(self.builddirectory, self.filename+self.ext), \
                        "-DCMAKE_BUILD_TYPE=RelWithDebInfo", \
                        "-DLOGLEVEL=INFO", \
                        "-DEXECUTABLE="+self.filename, \
                        os.path.join(self.tfdirectory,os.curdir)]
    for i in range(2):
      self.runcommand(command, dirname, exception=SimulationsErrorConfigure, logfilename="cmake."+repr(i)+".log")

  def build(self, force=False):

    dirname = os.path.join(self.builddirectory, "build")
    if force:
      command = ["make", "clean"]
      self.runcommand(command, dirname, exception=SimulationsErrorBuild)

    requiredinput = self.getrequiredinput(0, build=True)

    input_changed = False
    for filepath_k, filepath_v in requiredinput.items():
      try:
        try:
          fh = open(os.path.join(dirname, os.path.basename(filepath_v))).read()
        except UnicodeDecodeError:
          fh = gzip.open(os.path.join(dirname, os.path.basename(filepath_v)), 'rt', encoding='utf-8').read()
        checksum = hashlib.md5(fh.encode('utf-8')).hexdigest()
      except:
        checksum = None
      try:
        try:
          fh = open(filepath_k).read()
        except UnicodeDecodeError:
          fh = gzip.open(filepath_k, 'rt', encoding='utf-8').read()
        input_changed = input_changed or checksum != hashlib.md5(fh.encode('utf-8')).hexdigest()
      except IOError:
        self.log("WARNING: Unable to open %s"%(filepath_k))
        input_changed = True

    if input_changed or force:
      # file has changed or a recompilation is necessary
      for filepath_k, filepath_v in requiredinput.items():
        try:
          shutil.copy(filepath_k, os.path.join(dirname, os.path.basename(filepath_v)))
        except IOError:
          self.log("WARNING: required input (%s) not found, continuing anyway."%(filepath_k))

    env = copy.deepcopy(os.environ)
    try:
      env["PYTHONPATH"] = ":".join([dirname, env["PYTHONPATH"]])
    except KeyError:
      env["PYTHONPATH"] = dirname

    self.log("Building in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))
    command = ["make"]
    self.runcommand(command, dirname, exception=SimulationsErrorBuild, logfilename="make.log", verbose=True)

  def cleanbuild(self):
    try:
      shutil.rmtree(os.path.join(self.basedirectory, self.filename+self.ext+".build"))
    except OSError:
      pass

  def clean(self):
    Run.clean(self)
    self.cleanbuild()

  def numbercheckpoints(self, run=0):
    """Helper function to return the checkpoint directory and number of checkpoints."""
    
    dirname = os.path.join(self.rundirectory, "run_"+repr(run).zfill(len(repr(self.nruns))))

    depth = 0
    for root, dirnames, files in os.walk(dirname, topdown=True):
      depth = root.count(os.path.sep) - dirname.count(os.path.sep)
      if "checkpoint" in dirnames:
        # there is a checkpoint directory beneath the current directory
        dirnames[:] = ["checkpoint"]
      else:
        # there is no checkpoint directory beneath the current directory
        dirnames[:] = []

    return depth

  def checkpointrun(self, index=-1):

    for r in range(self.nruns):

      dirname = os.path.join(self.rundirectory, "run_"+repr(r).zfill(len(repr(self.nruns))))

      self.log("Checking for checkpoints in or beneath directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))

      try:
        threadlibspud.load_options(os.path.join(dirname, self.filename+self.ext))
        output_base_name = libspud.get_option("/io/output_base_name")
        threadlibspud.clear_options()
      except KeyError:
        self.log("  Unable to find output_base_name from base input file in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))
        continue

      depth = self.numbercheckpoints(run=r)

      checkpointdir = os.path.join("", *(["checkpoint"]*depth))
      basedir = os.path.join(dirname, checkpointdir)

      files = glob.glob1(basedir, "*"+self.ext)
      files = [f for f in files if output_base_name+"_checkpoint"*(depth+1) in f and f.split(self.ext)[0].split("_")[-1].isdigit()]
      files = sorted(files, key=lambda f: int(f.split(self.ext)[0].split("_")[-1]))

      try:
        basefile = files[index]
      except IndexError:
        self.log("  Requested checkpoint file not found in directory: %s"%(os.path.relpath(basedir, self.currentdirectory)))
        continue

      requiredinput = self.getrequiredcheckpointinput(basedir, basefile, r)
      commands = self.getcommands(basefile=basefile)

      dirname = os.path.join(basedir, "checkpoint")
      try:
        os.makedirs(dirname)
      except OSError:
        pass
        
      self.log("  Running from checkpoint in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))
      for filepath in requiredinput:
        shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
      self.writecheckpointoptions(basedir, basefile)
      
      for i, command in enumerate(commands):
        logf = '_'+self.filename+self.ext+'_checkpoint.'+repr(i)+'.log'
        self.runcommand(command, dirname, exception=SimulationsErrorRun, logfilename=logf)
      
      self.log("  Finished running from checkpoint in directory: %s"%(os.path.relpath(dirname, self.currentdirectory)))

  def getrequiredcheckpointinput(self, basedir, basefile, run):
    # get the standard required input for this simulation
    requiredinput = self.getrequiredinput(run)
    threadlibspud.load_options(os.path.join(basedir, basefile))
    for s in range(libspud.option_count("/system")):
      filename = libspud.get_option("/system["+repr(s)+"]/field/type/rank/initial_condition/file")
      # for checkpoints we assume that the checkpointed file is the one we want so 
      # we clean the requiredinput list of any references to it as a value (if any exist)
      popkeys = []
      for inputpath_k, inputpath_v in requiredinput.items():
        if filename == os.path.basename(inputpath_v): popkeys.append(inputpath_k)
      for popkey in popkeys: requiredinput.pop(popkey)
      requiredinput[os.path.join(basedir, filename)] = filename
    threadlibspud.clear_options()
    return requiredinput

  def getcommands(self, basefile=None):
    if basefile is None: basefile = self.filename+self.ext
    nprocs = self.getnprocs()
    valgrind_opts = self.optionsdict["valgrind"]
    mpi_opts = []
    if "mpi" in self.optionsdict: mpi_opts = [opt.strip() for opt in self.optionsdict["mpi"]]
    commands = [[]]
    if nprocs > 1:
      commands[0] += ["mpiexec"]+mpi_opts+["-np", repr(nprocs)]
    if valgrind_opts is not None:
      commands[0] += ["valgrind"]+valgrind_opts
    commands[0] += [os.path.join(self.builddirectory, "build", self.filename), "-vINFO", "-l", basefile]
    return commands

  def getnprocs(self):
    # take the base number of processes and scale it with all the values requested for the current parameters
    nprocs = reduce(operator.mul, list(self.optionsdict["procscales"].values()), self.optionsdict["nprocs"])
    return nprocs

  def getbuilddirectory(self):
    '''Return the path to the build directory for this simulation.'''

    # start at the base directory
    builddirectory = os.path.join(self.basedirectory, self.filename+self.ext+".build")

    # decide the order (preserve what exists if ordered dictionary or use alphabetical)
    if isinstance(self.optionsdict["values"], collections.OrderedDict):
      sorteditems = list(self.optionsdict["values"].items())
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
    for simulationpath, simoptionsdict in self.globaloptionsdict.items():
      # set up a list to contain the options dictionaries for each parameter set
      slicedoptionslist = []
      self.sliceoptionsdict(slicedoptionslist, simoptionsdict) # slice the global options dictionary parameter values
      for optionsdict in slicedoptionslist:
        # take a copy of the simulation options dictionary and set the 
        # values to the subset we've just extracted from the parameter sweep
        localsimoptionsdict = copy.deepcopy(simoptionsdict)
        localsimoptionsdict["values"]        = optionsdict["values"]
        localsimoptionsdict["procscales"]    = optionsdict["procscales"]
        localsimoptionsdict["value_indices"] = optionsdict["value_indices"]
        localsimoptionsdict["value_lengths"] = optionsdict["value_lengths"]
        # append the simulation to the list
        self.simulations.append(simoptionsdict["type"](simulationpath, localsimoptionsdict, \
                                                       currentdirectory, tfdirectory, dirname, \
                                                       logprefix=self.logprefix))

    self.runs = {}
    self.builds = {}
    self.collate()

    self.threadruns = []
    self.threadbuilds = []

    self.tests = [Test(name, code) for name, code in tests.items()]

  def sliceoptionsdict(self, slicedoptionslist, diroptionsdict, optionsdict=None, i=0):
    '''Recursively slice the global options dictionary into a list of simulation optionsdictionaries 
       that contain a single parameter set from the full parameter sweep.'''

    # if this is the first level then set up an ordered dictionary 
    # (this gets repeatedly reused but we take deep copies into the list to ensure there's no overwriting)
    if optionsdict is None: optionsdict = {"values"        : collections.OrderedDict(), \
                                           "procscales"    : collections.OrderedDict(),
                                           "value_indices" : collections.OrderedDict(),
                                           "value_lengths" : collections.OrderedDict() }

    assert len(diroptionsdict["values"]) == len(diroptionsdict["procscales"])

    # if we haven't reached the end yet...
    if i < len(diroptionsdict["values"]):

      # fetch the next sorted valitem
      if isinstance(diroptionsdict["values"], collections.OrderedDict):
        valitem = list(diroptionsdict["values"].items())[i]
      else:
        valitem = sorted(diroptionsdict["values"].items())[i]

      # fetch the next sorted procitem
      if isinstance(diroptionsdict["procscales"], collections.OrderedDict):
        procitem = list(diroptionsdict["procscales"].items())[i]
      else:
        procitem = sorted(diroptionsdict["procscales"].items())[i]

      assert procitem[0] == valitem[0]    # val and procitem [0] should contain the name of the same parameter
      # each entry of procitem[1] corresponds to the process scale for each entry of
      # valitem[1] (the values of the parameter) hence they should have the same length
      if len(procitem[1]) != len(valitem[1]):
        self.log("ERROR: number of process scales does not equal the number of values:")
        self.log("parameter %(parameter_name)s"%{"parameter_name": valitem[0]})
        raise SimulationsErrorInitialization

      optionsdict["value_lengths"][valitem[0]] = len(valitem[1])
      # loop over the values provided in the list associated with valitem
      for j in range(len(valitem[1])):
        # remember the index
        optionsdict["value_indices"][valitem[0]] = j
        # set the optionsdict to one of the values
        optionsdict["values"][valitem[0]] = valitem[1][j]
        # set the optionsdict to one of the procscales
        optionsdict["procscales"][procitem[0]] = procitem[1][j]
        # and recurse to the next level of the global options dictionary
        self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1)

    # we've reached the lowest level of the diroptionsdict["values"]
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
        # get the previous definition
        oldsim = self.runs[self.simulations[s].rundirectory]
        # pop the new definition
        tmpsim = self.simulations.pop(s)
        # join the variables of the old and new simulations
        oldsimvar = set([(var.name, var.code) for var in oldsim['simulation'].variables])
        oldsimvar = oldsimvar.union(set([(var.name, var.code) for var in tmpsim.variables]))
        oldsim['simulation'].variables = [Variable(var[0], var[1]) for var in oldsimvar]
        # join the run when tests of the old and new simulations
        for k,v in oldsim['simulation'].optionsdict['run_when'].items():
          if k in ["never"]:
            oldsim['simulation'].optionsdict['run_when'][k] = v and tmpsim.optionsdict['run_when'][k]
          else:
            oldsim['simulation'].optionsdict['run_when'][k] = v or tmpsim.optionsdict['run_when'][k]
        # delete the duplicate definition
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
    for run in self.runs.values():
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
    for s in range(len(simulation.dependencies)):
      # if the dependency rundirectory is already in the runs dictionary then we need to delete it and check its (d)level
      if simulation.dependencies[s].rundirectory in self.runs:
        # get the previous definition
        oldsim = self.runs[simulation.dependencies[s].rundirectory]
        # pop the new definition
        tmpsim = simulation.dependencies.pop(s)
        # add dependents of the new definition to the old definition
        oldsim['simulation'].dependents += tmpsim.dependents
        # join the variables of the old and new simulations
        oldsimvar = set([(var.name, var.code) for var in oldsim['simulation'].variables])
        oldsimvar = oldsimvar.union(set([(var.name, var.code) for var in tmpsim.variables]))
        oldsim['simulation'].variables = [Variable(var[0], var[1]) for var in oldsimvar]
        # join the run when tests of the old and new simulations
        for k,v in oldsim['simulation'].optionsdict['run_when'].items():
          if k in ["never"]:
            oldsim['simulation'].optionsdict['run_when'][k] = v and tmpsim.optionsdict['run_when'][k]
          else:
            oldsim['simulation'].optionsdict['run_when'][k] = v or tmpsim.optionsdict['run_when'][k]
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

  def configure(self, level=None, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator(self.simulationselector(self.builds, level=level, dlevel=dlevel, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadconfigure, args=[queue], kwargs={'force':force}), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def threadbuild(self, queue, force=False):
    error = None
    for simulation in self.threadbuilds: 
      try:
        simulation.build(force=force)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def build(self, level=None, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator(self.simulationselector(self.builds, level=level, dlevel=dlevel, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadbuild, args=[queue], kwargs={'force':force}), queue]) 
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def threadrun(self, queue, force=False, rundependencies=True):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.run(force=force, rundependencies=rundependencies)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def run(self, level=None, dlevel=0, types=None, force=False):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadrun, args=[queue], kwargs={'force':force, 'rundependencies':level==None}), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join()
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

    dlevel += 1
    # we request level=None here to make sure we recurse to all dlevels
    if len(self.simulationselector(self.runs, level=None, dlevel=dlevel, types=types)) > 0:
      self.writeoptions(level=level, dlevel=dlevel, types=types)
      self.configure(level=level, dlevel=dlevel, types=types, force=force)
      self.build(level=level, dlevel=dlevel, types=types, force=force)
      self.run(level=level, dlevel=dlevel, types=types, force=force)

  def threadcheckpointrun(self, queue, index=-1):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.checkpointrun(index=index)
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def checkpointrun(self, index=-1, level=None, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadcheckpointrun, args=[queue], kwargs={'index':index}), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def listinput(self, level=None, dlevel=0, types=None):
    runs = self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)[::-1]
    inputlist = set()
    for simulation in runs: inputlist.update(set(simulation.listinput()))
    return list(inputlist)

  def writeoptions(self, level=None, dlevel=0, types=None):
    runs = self.simulationselector(self.runs, level=level, dlevel=dlevel, types=types)[::-1]
    for simulation in runs: simulation.writeoptions()  # not thread safe so don't thread it!

  def threadcleanrun(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.cleanrun()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def cleanrun(self, level=None, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadcleanrun, args=[queue]), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def threadcleanbuild(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.cleanbuild()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def cleanbuild(self, level=None, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadcleanbuild, args=[queue]), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def threadclean(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.clean()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def clean(self, level=None, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadclean, args=[queue]), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def threadcheckpointclean(self, queue):
    error = None
    for simulation in self.threadruns: 
      try:
        simulation.checkpointclean()
      except:
        ex_type, ex_value, tb = sys.exc_info()
        error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
    queue.put(error)

  def checkpointclean(self, level=None, types=None):
    threadlist=[]
    self.threadruns = ThreadIterator(self.simulationselector(self.runs, level=level, types=types))
    for i in range(self.nthreads):
      queue = Queue()
      threadlist.append([threading.Thread(target=self.threadcheckpointclean, args=[queue]), queue])
      threadlist[-1][0].start()
    for t in threadlist:
      # wait until all threads finish
      t[0].join() 
    for t in threadlist:
      error = t[1].get()
      if error is not None:
        ex_type, ex_value, tb_str = error
        message = '%s (in thread)%s%s' % (str(ex_value), os.linesep, tb_str)
        raise Exception(message) 

  def evaluatevariables(self, level=None, dlevel=None, types=None):
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

      for name, value in simvarsdict.items():
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
            def dependentwalk(simulation):
              paths = []
              for dependent in simulation.dependents:
                paths += dependentwalk(dependent)
              if len(simulation.dependents)==0 or simulation.path in self.globaloptionsdict: 
                paths += [simulation.path]
              return paths
            paths = dependentwalk(simulation)
            values = collections.OrderedDict()
            for path in paths:
              values.update(self.globaloptionsdict[path]["values"])
            varsdict[name] = NestedList(values)
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
      try:
        status = test.run(varsdict)
      except TestOrVariableException as exc:
        self.log("ERROR: Evaluation of %s test raised an exception."%(str(test.name)))
        self.log(exc.message)
        self.log("failure.")
        teststatus.append('F')
      else:
        if status:
          self.log("success.")
          teststatus.append('P')
        else:
          self.log("failure.")
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
    for line in string.splitlines():
      if self.logprefix is not None: line = self.logprefix+": "+line
      sys.stdout.write(line+os.linesep)
    sys.stdout.flush()

  def simulationselector(self, simdict, level=None, dlevel=None, types=None):
    return [value['simulation'] \
            for value in sorted(list(simdict.values()), key=lambda value: value['level']) \
            if (level is None or value['level'] == level) and \
               (dlevel is None or value['dlevel'] == dlevel) and \
               (types is None or isinstance(value['simulation'], tuple(types)))]

#########################################################################################

class SimulationHarnessBatch(SimulationBatch):
  '''A derived SimulationBatch that takes in a list of spud based harnessfiles and sets up
     subgroups of SimulationBatches based on them.'''

  def __init__(self, harnessfiles, filename, currentdirectory, tfdirectory, nthreads=1, parameters={}, extraoptions={}):
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
    # any parameter values we want to overload
    self.parameters = parameters

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
      for d in range(libspud.option_count("/simulations/simulation")):
        simulation_optionpath = "/simulations/simulation["+repr(d)+"]"
        harnessfileoptionsdict.update(self.getoptions(simulation_optionpath, dirname, extraoptions=extraoptions))

      # loop over any runs
      for d in range(libspud.option_count("/simulations/run")):
        run_optionpath = "/simulations/run["+repr(d)+"]"
        harnessfileoptionsdict.update(self.getoptions(run_optionpath, dirname, run=True, extraoptions=extraoptions))

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
       sys.stdout.write(os.linesep)
       sys.stdout.write("Summary of test problems with failures:"+os.linesep)
       for group in failinggroups:
          sys.stdout.write(group[0]+': '+group[1]+os.linesep)
       sys.stdout.write(os.linesep)
       sys.stdout.flush()
    
    # report the overall number of passes and failures
    if passcount + failcount > 0:
      self.log("Passes:   %d" % passcount)
      self.log("Failures: %d" % failcount)
    
    # exit with an error if any failed
    if varerror: raise SimulationsErrorVariable
    if failcount > 0: raise SimulationsErrorTest

  def getrequiredfiles(self, optionpath, dirname=None, build=False):
     '''Get a list of the required files listed under the optionpath.  Run python scripts in dirname if supplied.'''

     # regex to split up filename lists
     # accept comma, semicolon, space or newline delimiters
     # NODE                     EXPLANATION
     # --------------------------------------------------------------------------------
     #   (?:                      group, but do not capture (1 or more times
     #                            (matching the most amount possible)):
     # --------------------------------------------------------------------------------
     #     [^,; \n]                 any character except: ',', ';', ' ',
     #                              '\n' (newline)
     # --------------------------------------------------------------------------------
     #   )+                       end of grouping
     # - from http://rick.measham.id.au/paste/explain.pl
     r = re.compile(r'(?:[^,; '+os.linesep+'])+')

     required_files = []
     required_files_python = []
     # loop over the filenames
     for f in range(libspud.option_count(optionpath+"/filenames")):
       filename_optionpath = optionpath+"/filenames["+repr(f)+"]"
       if build and not libspud.have_option(filename_optionpath+"/required_at_build"): continue
       # try getting the files as a string
       try:
         filenames = libspud.get_option(filename_optionpath+"/string")
         required_files += r.findall(filenames)
       # otherwise run a python string to get the filename(s)
       except libspud.SpudKeyError:
         # can't do anything with these yet as we need the simulation
         # options dict as an environment
         filename_name = libspud.get_option(filename_optionpath+"/name")
         filename_code = libspud.get_option(filename_optionpath+"/python")
         required_files_python.append((filename_name, filename_code))
         
     # if dirname has been supplied then make the filenames relative to it
     if dirname is not None:
       for f in range(len(required_files)):
         filename = required_files[f]
         if os.path.isabs(filename):
           required_files[f] = os.path.normpath(filename)
         else:
           required_files[f] = os.path.normpath(os.path.join(dirname, filename))

     # return the required files
     return required_files, required_files_python

  def getcommands(self, optionpath):
     commands = []
     for c in range(libspud.option_count(optionpath+"/command")):
       command_optionpath = optionpath+"/command["+repr(c)+"]"
       commands.append(libspud.get_option(command_optionpath).split())
     return commands

  def getvariables(self, optionpath):
     variables  = collections.OrderedDict()
     for v in range(libspud.option_count(optionpath+"/variable")):
       variable_optionpath = optionpath+"/variable["+repr(v)+"]"
       variable_name = libspud.get_option(variable_optionpath+"/name")
       variables[variable_name] = libspud.get_option(variable_optionpath)
     return variables

  def gettests(self, optionpath):
     tests  = collections.OrderedDict()
     for t in range(libspud.option_count(optionpath+"/test")):
       test_optionpath = optionpath+"/test["+repr(t)+"]"
       test_name = libspud.get_option(test_optionpath+"/name")
       tests[test_name] = libspud.get_option(test_optionpath)
     return tests

  def getparameters(self, optionpath, parent_parameters=None):
     '''Fetch parameters from /parameters child of the given optionpath.  
        Test that the parameters are available in the parent if parent_parameters present.'''

     # regular expression to split up values string into a list
     # the list may be comma(,), semicolon(;), space ( ) or newline (\n) delimited
     # if the list is of bracketed items (e.g. tuples) then any delimiters within 
     # the brackets are preserved
     # brackets may be (), [] or {}
     #NODE                     EXPLANATION
     #--------------------------------------------------------------------------------
     #  (?:                      group, but do not capture (1 or more times
     #                           (matching the most amount possible)):
     #--------------------------------------------------------------------------------
     #    [^,; \n([{]              any character except: ',', ';', ' ',
     #                             '\n' (newline), '(', '[', '{'
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
     r = re.compile(r'(?:[^,; '+os.linesep+'([{]|\([^)]*\)|\[[^]]*\]|\{[^}]*\})+')
    
     # set up ordered dictionaries to remember the order in which parameters are specified
     parameter_values    = collections.OrderedDict()
     parameter_updates   = collections.OrderedDict()
     parameter_builds    = collections.OrderedDict()
     parameter_procscales = collections.OrderedDict()
     # loop over all the parameters
     for p in range(libspud.option_count(optionpath+"/parameter")):
       parameter_optionpath = optionpath+"/parameter["+repr(p)+"]"
       parameter_name = libspud.get_option(parameter_optionpath+"/name")
       
       # throw an error if the parameter is not in the parent list (if present)
       if parent_parameters is not None:
         if parameter_name not in parent_parameters:
           self.log("ERROR: dependency parameter not in parent simulation:")
           self.log("parameter %(parameter_name)s (%(parameter_optionpath)s)"%\
                 {"parameter_name": parameter_name, "parameter_optionpath": parameter_optionpath})
           
           raise SimulationsErrorInitialization
       
       # try to get a list of the parameter values (not necessarily specified for child parameters so except failure)
       if libspud.have_option(parameter_optionpath+"/values"):
         if parameter_name in self.parameters:
           parameter_values[parameter_name] = self.parameters[parameter_name]
         else:
           parameter_values[parameter_name]  = r.findall(libspud.get_option(parameter_optionpath+"/values"))
         # try to get a list of the process scales (not necessarily specified so except failure)
         try:
           parameter_procscales[parameter_name] = libspud.get_option(parameter_optionpath+"/process_scale")
         except libspud.SpudKeyError:
           parameter_procscales[parameter_name] = [1 for v in range(len(parameter_values[parameter_name]))]
       else:
         # insert nothing if we haven't specified any values
         # try to get a list of the process scales (not necessarily specified so except failure)
         try:
           parameter_procscales[parameter_name] = libspud.get_option(parameter_optionpath+"/process_scale")
         except libspud.SpudKeyError:
           parameter_procscales[parameter_name] = None # haven't specified how we want processes scaled (defaults to 1 eventually)

       # try to get a list of the parameter update instructions (not necessarily specified for parent parameters so except failure)
       try:
         parameter_updates[parameter_name] = libspud.get_option(parameter_optionpath+"/update")
         parameter_builds[parameter_name]  = libspud.have_option(parameter_optionpath+"/update/single_build")
       except libspud.SpudKeyError:
         parameter_updates[parameter_name] = None # insert None if we haven't said how to update
         parameter_builds[parameter_name]  = True # not updating so we don't need multiple builds

     # return the OrderedDict objects for the values and updates
     return parameter_values, parameter_updates, parameter_builds, parameter_procscales

  def getdependencies(self, optionpath, dirname, parent_parameters):
     dependencies_options = {}

     for d in range(libspud.option_count(optionpath+"/simulation")):
        simulation_optionpath = optionpath+"/simulation["+repr(d)+"]"
        dependencies_options.update(self.getoptions(simulation_optionpath, dirname, parent_parameters))
     
     for d in range(libspud.option_count(optionpath+"/run")):
        run_optionpath = optionpath+"/run["+repr(d)+"]"
        dependencies_options.update(self.getoptions(run_optionpath, dirname, parent_parameters, run=True))
     
     return dependencies_options

  def getoptions(self, optionpath, dirname, parent_parameters=None, run=False, extraoptions={}):
     options = {}
     name = libspud.get_option(optionpath+"/name")

     # name of options file and normalize path to it relative to the path to the harnessfile itself
     input_file = libspud.get_option(optionpath+"/input_file")
     if os.path.isabs(input_file):
       path = os.path.normpath(input_file)
     else:
       path = os.path.normpath(os.path.join(dirname, input_file))
     
     try:
       run_when = libspud.get_option(optionpath+"/run_when/name")
     except libspud.SpudKeyError:
       run_when = "input_changed_or_output_missing"

     # fetch required input and output for this simulation (used to copy files into rundirectory 
     # and establish if run has successfully completed
     required_input, required_input_python  = \
                        self.getrequiredfiles(optionpath+"/required_input", dirname)
     required_build_input, required_build_input_python  = \
                        self.getrequiredfiles(optionpath+"/required_input", dirname, build=True)
     required_output, required_output_python = \
                        self.getrequiredfiles(optionpath+"/required_output")

     # get the dictionary of the parameters we are sweeping over in this simulation and the number of times each run should
     # be performed (e.g. for benchmarking purposes)
     parameter_values, parameter_updates, \
                       parameter_builds,  \
                       parameter_procscales = self.getparameters(optionpath+"/parameter_sweep", \
                                                                parent_parameters=parent_parameters)
     try:
       nruns = libspud.get_option(optionpath+"/parameter_sweep/number_runs")
     except libspud.SpudKeyError:
       nruns = 1
     
     variables = self.getvariables(optionpath+"/variables")

     options[path] = {}
     options[path].update(extraoptions)
     options[path]["name"]      = name
     options[path]["run_when"] = {                                                         \
                                             "input_changed"  : run_when.find("input_changed")  != -1, \
                                             "output_missing" : run_when.find("output_missing") != -1, \
                                             "always"         : run_when == "always",                  \
                                             "never"          : run_when == "never",                   \
                                            }
     if run:
       options[path]["type"]    = Run
       options[path]["spudfile"] = libspud.have_option(optionpath+"/input_file/spud_file")
       options[path]["run"] = self.getcommands(optionpath+"/commands")
     else:
       options[path]["type"]    = Simulation
       try:
         options[path]["nprocs"] = libspud.get_option(optionpath+"/number_processes")
       except libspud.SpudKeyError:
         options[path]["nprocs"] = 1
       try:
         options[path]["valgrind"] = libspud.get_option(optionpath+"/valgrind_options").split()
       except libspud.SpudKeyError:
         options[path]["valgrind"] = None
       # get the parameters for any checkpoint pickups
       checkpoint_values, checkpoint_updates, \
                          checkpoint_builds,  \
                          checkpoint_procscales = self.getparameters(optionpath+"/checkpoint")
       # don't do anything with checkpoint_builds or checkpoint_procscales
       # but make sure we only have single values and take them out of list form
       for k, v in checkpoint_values.items():
         if len(v) > 1:
           self.log("ERROR: cannot have more than one value for checkpoint parameters.")
           self.log("parameter: %(parameter_name)s"%{"parameter_name": k})
           raise SimulationsErrorInitialization
         checkpoint_values[k] = v[0]
       options[path]["checkpoint_values"]  = checkpoint_values
       options[path]["checkpoint_updates"] = checkpoint_updates
     options[path]["input"]              = required_input
     options[path]["input_python"]       = required_input_python
     options[path]["build_input"]        = required_build_input
     options[path]["build_input_python"] = required_build_input_python
     options[path]["output"]             = required_output
     options[path]["output_python"]      = required_output_python
     options[path]["values"]             = parameter_values
     options[path]["updates"]       = parameter_updates
     options[path]["builds"]        = parameter_builds
     options[path]["procscales"]    = parameter_procscales # not needed for runs but easier to keep in as dummy
     options[path]["nruns"] = nruns
     options[path]["variables"] = variables
     # if the simulation has dependencies then recursively get their options too
     if libspud.have_option(optionpath+"/dependencies"):
       options[path]["dependencies"] = \
                           self.getdependencies(optionpath+"/dependencies", dirname, \
                                                list(parameter_updates.keys()))

     return options

