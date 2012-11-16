import os
import subprocess
import sys
import hashlib
import shutil
import threading
import copy

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

class SimulationBatch:
  def __init__(self, name, ext, scriptdirectory, basebucketdirectory, globaloptionsdict, simulationtype=None, nthreads=1):
    self.globaloptionsdict = globaloptionsdict
    self.nthreads = nthreads
    self.simulations = []
    for dirname, diroptionsdict in self.globaloptionsdict.iteritems():
      slicedoptionslist = []
      self.sliceoptionsdict(slicedoptionslist, diroptionsdict["options"])
      for optionsdict in slicedoptionslist:
        if simulationtype==None:
          self.simulations.append(Simulation(name, ext, scriptdirectory, dirname, basebucketdirectory, optionsdict))
        else:
          self.simulations.append(simulationtype(name, ext, scriptdirectory, dirname, basebucketdirectory, optionsdict))
    self.runs = {}
    self.collateruns()
    self.builds = {}
    self.collatebuilds() # has to be called after collateruns
    self.threadruns = []
    self.threadbuilds = []
  
  def sliceoptionsdict(self, slicedoptionslist, diroptionsdict, optionsdict={}, i=0):
    if i==0: optionsdict = {}
    if i < len(diroptionsdict): 
      item = sorted(diroptionsdict.items())[i]
      for strvalue in item[1]:
        optionsdict[item[0]] = strvalue
        self.sliceoptionsdict(slicedoptionslist, diroptionsdict, optionsdict=optionsdict, i=i+1)
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

  def writeoptions(self, level=-1):
    runs = [value['simulation'] for value in sorted(self.runs.values(), key=lambda value: value['level']) if value['level'] >= level]
    for simulation in runs: simulation.writeoptions()  # not thread safe so don't thread it!

class Simulation:
  def __init__(self, name, ext, scriptdirectory, basedirectory, basebucketdirectory, optionsdict, nruns=1, dependents=[]):
    self.optionsdict = optionsdict
    self.name = name
    self.ext = ext
    self.scriptdirectory = scriptdirectory
    self.basebucketdirectory = basebucketdirectory
    self.basedirectory = os.path.join(self.scriptdirectory, basedirectory)  # directory in which the options file is stored
    self.rundirectory = self.basedirectory
    for directory in [item[0]+"_"+item[1] for item in sorted(optionsdict.items())]: self.rundirectory = os.path.join(self.rundirectory, directory)
    self.builddirectory = self.rundirectory
    self.nruns = nruns
    self.dependencies = []
    self.dependents = dependents
    self.lock=threading.Lock()

  def writeoptions(self):
    print "ERROR: writeoptions should be overloaded."
    sys.exit(1)

  def updateoptions(self):
    print "ERROR: updateoptions should be overloaded."
    sys.exit(1)

  def configure(self, force=False):
    
    self.lock.acquire()

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

    self.lock.release()

  def build(self, force=False):

    self.lock.acquire()

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

    self.lock.release()

  def run(self, force=False):
    
    self.lock.acquire()
    
    for dependency in self.dependencies: dependency.run(force=force)
    
    print "checking in directory ", self.rundirectory

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
        print "  running in directory ", self.rundirectory
        # file has changed or a recompilation is necessary
        for filepath in requiredinput:
          shutil.copy(filepath, os.path.join(dirname, os.path.basename(filepath)))
        
        for command in commands:
          p = subprocess.Popen(command, cwd=dirname)
          retcode = p.wait()
          if retcode!=0:
            print "ERROR Command ", command, " returned ", retcode, " in directory: ", dirname
            #sys.exit(1)
        
        print "  finished in directory ", self.rundirectory

    self.lock.release()
      
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

  def getcommands(self):
    commands = [[os.path.join(self.builddirectory, "build", self.name), "-vINFO", "-l", self.name+self.ext]]
    return commands

