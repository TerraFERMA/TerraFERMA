import os
import subprocess
import sys
import hashlib
import shutil
import threading
import copy
import string
import glob

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

class Simulation:
  def __init__(self, name, ext, scriptdirectory, basedirectory, basebucketdirectory, \
               optionsdict, nruns=1, dependents=[], checkpointdict={}):
    self.optionsdict = optionsdict
    self.checkpointdict = checkpointdict
    self.name = name
    self.ext = ext
    self.scriptdirectory = scriptdirectory
    self.basebucketdirectory = basebucketdirectory
    self.basedirectory = os.path.join(self.scriptdirectory, basedirectory)  # directory in which the options file is stored
    self.rundirectory = rundirectorystr(self.basedirectory, optionsdict)
    self.builddirectory = self.rundirectory
    self.nruns = nruns
    self.dependencies = []
    self.dependents = dependents
    self.lock=threading.Lock()

  def writeoptions(self):
    print "ERROR: writeoptions should be overloaded."
    sys.exit(1)

  def writecheckpointoptions(self):
    print "ERROR: writecheckpointoptions should be overloaded."
    sys.exit(1)

  def updateoptions(self, optionsdict):
    print "ERROR: updateoptions should be overloaded."
    sys.exit(1)

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
                          "-DOPTIONSFILE="+os.path.join(self.builddirectory, self.name+self.ext), \
                          "-DCMAKE_BUILD_TYPE=RelWithDebInfo", \
                          "-DLOGLEVEL=INFO", \
                          "-DEXECUTABLE="+self.name, \
                          os.path.join(self.basebucketdirectory,os.curdir)],
                          cwd=dirname)
    retcode = p.wait()
    if retcode!=0:
      print "ERROR cmake returned ", retcode, " in directory: ", dirname
      sys.exit(1)

  def build(self, force=False):

    dirname = os.path.join(self.builddirectory, "build")
    if force:
      p = subprocess.Popen(["make", "clean"], cwd=dirname)
      retcode = p.wait()
      if retcode!=0:
        print "ERROR make clean returned ", retcode, " in directory: ", dirname
        sys.exit(1)

    p = subprocess.Popen(["make"], cwd=dirname)
    retcode = p.wait()
    if retcode!=0:
      print "ERROR make returned ", retcode, " in directory: ", dirname
      sys.exit(1)

  def checkpointrun(self):

    if len(self.checkpointdict)==0: return None
    
    for r in range(self.nruns):
      dirname = os.path.join(self.rundirectory, "run_"+`r`)

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

      print "checking in directory ", os.path.relpath(basedir, self.scriptdirectory)

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
        
      print "  running in directory ", os.path.relpath(dirname, self.scriptdirectory)
      for filepath in requiredinput:
        shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
      self.writecheckpointoptions(basedir, basefile)
      
      for command in commands:
        p = subprocess.Popen(command, cwd=dirname)
        retcode = p.wait()
        if retcode!=0:
          print "ERROR Command ", command, " returned ", retcode, " in directory: ", dirname
          #sys.exit(1)
      
      print "  finished in directory ", os.path.relpath(dirname, self.scriptdirectory)

  def run(self, force=False):
    
    self.lock.acquire()
    
    for dependency in self.dependencies: dependency.run(force=force)
    
    print "checking in directory ", os.path.relpath(self.rundirectory, self.scriptdirectory)

    requiredinput = self.getrequiredinput()
    requiredoutput = self.getrequiredoutput()
    commands = self.getcommands()

    for r in range(self.nruns):
      dirname = os.path.join(self.rundirectory, "run_"+`r`)
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
        input_changed = input_changed or checksum != hashlib.md5(open(filepath).read()).hexdigest()

      output_missing = False
      for filepath in requiredoutput:
        try:
          output_file = open(os.path.join(dirname, filepath))
          output_file.close()
        except IOError:
          output_missing = True
      if len(requiredoutput)==0: output_missing=True # don't know what output is needed so we have to force
        
      if output_missing or input_changed or force:
        print "  running in directory ", os.path.relpath(dirname, self.scriptdirectory)
        # file has changed or a recompilation is necessary
        for filepath in requiredoutput:
          try:
            os.remove(os.path.join(dirname, filepath))
          except OSError:
            pass
        for filepath in requiredinput:
          shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
        
        for command in commands:
          p = subprocess.Popen(command, cwd=dirname)
          retcode = p.wait()
          if retcode!=0:
            print "ERROR Command ", command, " returned ", retcode, " in directory: ", dirname
            #sys.exit(1)
        
        print "  finished in directory ", os.path.relpath(dirname, self.scriptdirectory)

    self.lock.release()
      
  def postprocess(self):
    return None

  def extract(self, options):
    suboptionsdict = {}
    for option in options: suboptionsdict[option] = self.optionsdict[option]
    return suboptionsdict

  def getrequiredoutput(self):
    requiredoutput = []
    return requiredoutput
 
  def getrequiredinput(self):
    requiredinput = [os.path.join(self.rundirectory, self.name+self.ext)]
    return requiredinput

  def getcheckpointcommands(self, basefile):
    commands = [[os.path.join(self.builddirectory, "build", self.name), "-vINFO", "-l", basefile]]
    return commands

  def getcommands(self):
    commands = [[os.path.join(self.builddirectory, "build", self.name), "-vINFO", "-l", self.name+self.ext]]
    return commands

class SimulationBatch:
  def __init__(self, name, ext, scriptdirectory, basebucketdirectory, globaloptionsdict, \
               simulationtype=Simulation, nthreads=1):
    self.globaloptionsdict = globaloptionsdict
    self.scriptdirectory = scriptdirectory
    self.nthreads = nthreads
    self.simulations = []
    for dirname, diroptionsdict in self.globaloptionsdict.iteritems():
      slicedoptionslist = []
      self.sliceoptionsdict(slicedoptionslist, diroptionsdict["options"])
      checkpointdict = {}
      if "checkpoint" in diroptionsdict: checkpointdict = diroptionsdict["checkpoint"]
      for optionsdict in slicedoptionslist:
        self.simulations.append(simulationtype(name, ext, scriptdirectory, dirname, basebucketdirectory, \
                                               optionsdict, checkpointdict=checkpointdict))
    self.runs = {}
    self.collateruns()
    self.builds = {}
    self.collatebuilds() # has to be called after collateruns
    self.threadruns = []
    self.threadbuilds = []
  
  def sliceoptionsdict(self, slicedoptionslist, diroptionsdict, optionsdict={}, i=0, orgoptiondict=None):
    if i==0: optionsdict = {}
    if i < len(diroptionsdict):
      item = sorted(diroptionsdict.items())[i]
      if orgoptiondict is not None:
        if item[0] in orgoptiondict:
          self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1, orgoptiondict=orgoptiondict)
        else:
          for strvalue in item[1]:
            optionsdict[item[0]] = strvalue
            self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1, orgoptiondict=orgoptiondict)
      else:
        for strvalue in item[1]:
          optionsdict[item[0]] = strvalue
          self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1, orgoptiondict=orgoptiondict)
    else:
      if orgoptiondict is not None:
        assert(len(orgoptiondict)==1)
        item = orgoptiondict.items()[0]
        if item[0] in diroptionsdict:
          orgoptionslist = []
          for strvalue in item[1]:
            if strvalue in diroptionsdict[item[0]]:
              optionsdict[item[0]] = strvalue
              orgoptionslist.append(copy.deepcopy(optionsdict))
            else:
              orgoptionslist.append(None)
          slicedoptionslist.append(orgoptionslist)
      else:
        slicedoptionslist.append(copy.deepcopy(optionsdict))

  def collatebuilds(self):
    self.builds = {}
    for run in self.runs.itervalues():
      self.builds[run['simulation'].builddirectory] = {'level':run['level'], 'simulation':run['simulation']}

  def collateruns(self):
    self.runs = {}
    s = 0
    while s < len(self.simulations):
      if self.simulations[s].rundirectory in self.runs:
        tmpsim = self.simulations.pop(s)
        del tmpsim
      else:
        self.runs[self.simulations[s].rundirectory] = {'level':0, 'simulation':self.simulations[s]}
        self.collatedependencies(self.simulations[s])
        s += 1

  def collatedependencies(self, simulation, level=1):
    for s in range(len(simulation.dependencies)):
      if simulation.dependencies[s].rundirectory in self.runs:
        oldsim = self.runs[simulation.dependencies[s].rundirectory]
        tmpsim = simulation.dependencies.pop(s)
        oldsim['simulation'].dependents += tmpsim.dependents
        del tmpsim
        simulation.dependencies.insert(s, oldsim['simulation'])
        oldsim['level'] = max(level, oldsim['level'])
      else:
        self.runs[simulation.dependencies[s].rundirectory] = {'level':level, 'simulation':simulation.dependencies[s]}
        self.collatedependencies(simulation.dependencies[s], level=level+1)

  def threadconfigure(self, force=False):
    for simulation in self.threadbuilds: simulation.configure(force=force)

  def configure(self, level=-1, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator([value['simulation'] for value in sorted(self.builds.values(), key=lambda value: value['level']) if value['level'] >= level])
    for i in range(self.nthreads):
      threadlist.append(threading.Thread(target=self.threadconfigure, kwargs={'force':force}))
      threadlist[-1].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t.join()

  def threadbuild(self, force=False):
    for simulation in self.threadbuilds: simulation.build(force=force)

  def build(self, level=-1, force=False):
    threadlist=[]
    self.threadbuilds = ThreadIterator([value['simulation'] for value in sorted(self.builds.values(), key=lambda value: value['level']) if value['level'] >= level])
    for i in range(self.nthreads):
      threadlist.append(threading.Thread(target=self.threadbuild, kwargs={'force':force})) 
      threadlist[-1].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t.join()

  def threadrun(self, force=False):
    for simulation in self.threadruns: simulation.run(force=force)

  def run(self, level=-1, force=False):
    threadlist=[]
    self.threadruns = ThreadIterator([value['simulation'] for value in sorted(self.runs.values(), key=lambda value: value['level']) if value['level'] >= level])
    for i in range(self.nthreads):
      threadlist.append(threading.Thread(target=self.threadrun, kwargs={'force':force}))
      threadlist[-1].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t.join() 

  def threadcheckpointrun(self):
    for simulation in self.threadruns: simulation.checkpointrun()

  def checkpointrun(self, level=-1):
    threadlist=[]
    self.threadruns = ThreadIterator([value['simulation'] for value in sorted(self.runs.values(), key=lambda value: value['level']) if value['level'] >= level])
    for i in range(self.nthreads):
      threadlist.append(threading.Thread(target=self.threadcheckpointrun))
      threadlist[-1].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t.join() 

  def threadpostprocess(self):
    for simulation in self.threadruns: simulation.postprocess()

  def postprocess(self, level=-1):
    threadlist=[]
    self.threadruns = ThreadIterator([value['simulation'] for value in sorted(self.runs.values(), key=lambda value: value['level']) if value['level'] >= level])
    for i in range(self.nthreads):
      threadlist.append(threading.Thread(target=self.threadpostprocess))
      threadlist[-1].start()
    for t in threadlist:
      '''Wait until all threads finish'''
      t.join() 

  def writeoptions(self, level=-1):
    runs = [value['simulation'] for value in sorted(self.runs.values(), key=lambda value: value['level']) if value['level'] >= level]
    for simulation in runs: simulation.writeoptions()  # not thread safe so don't thread it!

def rundirectorystr(basedirectory, optionsdict):
    rundirectory = basedirectory
    for directory in [item[0]+"_"+item[1] for item in sorted(optionsdict.items())]: rundirectory = os.path.join(rundirectory, directory)
    return rundirectory

