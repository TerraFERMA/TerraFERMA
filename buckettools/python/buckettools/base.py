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

import sys
import re
import hashlib
import shutil
import subprocess
import os

# Functions used to provide comments in ufls
def comment(string):
  """Returns a ufl formatted comment based around the provided string."""
  return "# "+string+os.linesep

def declaration_comment(thing, type, name):
  """Comments the declaration of an object in ufl."""
  return comment(thing+" declaration for "+type+": "+name)

def produced_comment():
  """Comments the genesis of a file."""
  return comment("Produced by: "+" ".join(sys.argv))

# Functions used to provide single lines of ufl.
# These do not live in a class because they are needed by more than
# one unrelated class (i.e. are more general than the description
# of any one class).
def equal_ufl(symbol_a, symbol_b, suffix=""):
  """Returns a ufl string setting two variables equal."""
  return symbol_a+suffix+" = "+symbol_b+suffix

def coefficient_ufl(symbol, suffix=""):
  """Returns a ufl string declaring a coefficient on an element."""
  return symbol+suffix+" = Coefficient("+symbol+"_e)"+os.linesep

def forms_ufl(form_symbols):
  """Returns a ufl string declaring a set of forms with potentially non-default symbols."""
  ufl = "forms = ["
  for s in range(len(form_symbols)-1):
    ufl += form_symbols[s]+", "
  ufl += form_symbols[-1]
  ufl += "]"+os.linesep
  return ufl

def testfunction_ufl(symbol):
  """Returns a ufl string declaring a test function on an element."""
  return symbol+"_t = TestFunction("+symbol+"_e)"+os.linesep

def trialfunction_ufl(symbol):
  """Returns a ufl string declaring a trial function on an element."""
  return symbol+"_a = TrialFunction("+symbol+"_e)"+os.linesep

def ufl_reserved():
  """Returns an (incomplete) array of ufl symbols that are considered reserved."""
  return ['rhs', 'lhs', 'action', 'dx', 'ds', 'dS', 'inner', 'dot', 'grad', 'div', 'sym', 'derivative']

def uflsymbol_suffixes():
  """Returns an array of available function symbols in the ufl."""
  suffixes = function_uflsymbol_suffixes()
  suffixes += ["_e", "_e0", "_e1"]
  return suffixes

def function_uflsymbol_suffixes():
  """Returns an array of available function symbols in the ufl."""
  return ["", "_i", "_n"]

def form_symbols(namespace):
   try:
     headerfile = open(namespace+".h", 'r')
   except IOError:
     sys.stdout.write("ERROR: unable to open %s for reading."%(namespace+".h")+os.linesep)
     sys.stdout.write("       This file must be generated before the form symbols can be ascertained."+os.linesep)
     sys.stdout.flush()
     sys.exit(1)

   headertext = headerfile.read()
   return [line[23:].split()[0][:-1] for line in headertext.splitlines() if line.startswith("class CoefficientSpace_")]

def ffc(namespace, form_representation, quadrature_rule, quadrature_degree):
  uflfilename = namespace+".ufl"

  try:
    checksum = hashlib.md5(open(uflfilename, 'r').read()).hexdigest()
  except:
    checksum = None

  try:
    headerfile = open(namespace+".h", 'r')
  except IOError:
    headerfile = None

  if headerfile:
    rebuild = checksum != hashlib.md5(open(uflfilename+".temp", 'r').read()).hexdigest()

    headertext = headerfile.read()
    
    qdi  = headertext.find("quadrature_degree")
    qdin = headertext.find(os.linesep, qdi, -1)
    if (qdi != -1) and (qdin != -1):
      qd_old = headertext[qdi:qdin].split(" ")[-1]
      qd_new = `quadrature_degree` if quadrature_degree is not None else "-1"
      rebuild = rebuild or qd_new != qd_old
    else: # if we don't know what the old quadrature degree was we always want to (re)build
      rebuild = True

    qri  = headertext.find("quadrature_rule")
    qrin = headertext.find(os.linesep, qri, -1)
    if (qri != -1) and (qrin != -1):
      qr_old = headertext[qri:qrin].split(" ")[-1]
      qr_new = "'"+quadrature_rule+"'" 
      rebuild = rebuild or qr_new != qr_old
    else: # if we don't know what the old quadrature rule was we always want to (re)build
      rebuild = True
    
    fri  = headertext.find("representation")
    frin = headertext.find(os.linesep, fri, -1)
    if (fri != -1) and (frin != -1):
      fr_old = headertext[fri:frin].split(" ")[-1]
      fr_new = "'"+form_representation+"'" 
      rebuild = rebuild or fr_new != fr_old
    else: # if we don't know what the old representation was we always want to (re)build
      rebuild = True
    
  else:   # if we don't have a header file we always want to (re)build
    rebuild = True

  if rebuild:
    # files and/or quadrature degree have changed
    shutil.copy(uflfilename+".temp", uflfilename)
    command = ["ffc", "-l", "dolfin", "-O", "-r", form_representation]
    if quadrature_degree is not None:
      command += ["-fquadrature_degree="+`quadrature_degree`]
    command += ["-fsplit"]
    command += ["-q", quadrature_rule]
    command += [uflfilename]
    logfilename = uflfilename[:-4]+".ffc.log"
    sys.stdout.write("Running: "+" ".join(command) + " > %s"%(logfilename) + os.linesep); sys.stdout.flush()
    # start running the ffc command
    try:
      p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception as exc:
      sys.stdout.write("ERROR: ffc raised an exception before running on file %s"%(uflfilename)+os.linesep); sys.stdout.flush()
      sys.stdout.write(str(exc))
      sys.exit(1)
    # search the output for warnings
    warning = re.compile("WARNING:")
    warnings = False
    output = ""  # save the output in case we want to pipe it to stdout on failure
    with open(logfilename, 'w') as logfile: # but otherwise just send it to a logfile
      i = -1
      for line in iter(p.stdout.readline, ""):
        logfile.write(line)
        output += line
        if warning.search(line):            # search for warnings
          sys.stdout.write(line)
          i = line.find("WARNING:")
          warnings = True
        elif i > -1 and line.startswith(" "*(i+8)):  # and indented lines after warnings
          sys.stdout.write(line)
        else:
          i = -1
    # report warnings
    if warnings:
      sys.stdout.write(os.linesep+"WARNING: ffc warnings detected!"+os.linesep*2)
    # wait until the command finishes (probably not necessary due to log above)
    retcode = p.wait()
    # if we've received a non-zero exit code
    if retcode:
      # report the error
      sys.stdout.write("ERROR: while calling ffc on file %s"%(uflfilename)+os.linesep); sys.stdout.flush()
      # and output everything we know
      for line in output.splitlines():
        sys.stdout.write(" "*2+line+os.linesep); sys.stdout.flush()
      # then check the output for hints on debugging
      verbose = "--verbose" in output
      debug = uflfilename[:-4]+"_debug.py" in output
      # if we were told to rerun with verbose enabled, then do it
      if verbose:
        command.insert(1, "--verbose")
        logfilename = uflfilename[:-4]+".ffc-verbose.log"
        sys.stdout.write(os.linesep)
        sys.stdout.write("Re-running with --verbose: " + " ".join(command) + " > %s"%(logfilename) + os.linesep)
        sys.stdout.flush()
        try:
          p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception as exc:
          sys.stdout.write("ERROR: ffc raised an exception before running on file %s"%(uflfilename)+os.linesep); sys.stdout.flush()
          sys.stdout.write(str(exc))
          sys.exit(1)
        with open(logfilename, 'w') as logfile:
          for line in iter(p.stdout.readline, ""):
            logfile.write(line)
            sys.stdout.write(" "*4+line); sys.stdout.flush()
        retcode = p.wait()
      # if a debugging script has been generated then run it
      if debug:
        command = ["python", uflfilename[:-4]+"_debug.py"]
        logfilename = uflfilename[:-4]+".ffc-debug.log"
        sys.stdout.write(os.linesep)
        sys.stdout.write("Running debugging script: " + " ".join(command) + " > %s"%(logfilename) + os.linesep)
        sys.stdout.flush()
        try:
          p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception as exc:
          sys.stdout.write("ERROR: python raised an exception before running on file %s"%(uflfilename[:-4]+"_debug.py")+os.linesep); sys.stdout.flush()
          sys.stdout.write(str(exc))
          sys.exit(1)
        with open(logfilename, 'w') as logfile:
          for line in iter(p.stdout.readline, ""):
            logfile.write(line)
            sys.stdout.write(" "*4+line); sys.stdout.flush()
        retcode = p.wait()
      # finally, don't let the failed input keep the same name (this can cause spurious successful builds on later attempts)
      sys.stdout.write(os.linesep)
      sys.stdout.write("Moving failed input to: %s"%(uflfilename+".fail")+os.linesep)
      shutil.move(uflfilename, uflfilename+".fail")
      sys.stdout.write(os.linesep)
      sys.exit(1)
    else:
      # we've succeeded at build (yay!)
      # on the off-chance we have previously tried to build but failed, let's clean up that output
      cleanupfilenames = [uflfilename+".fail", \
                          uflfilename[:-4]+"_debug.py", \
                          uflfilename[:-4]+".ffc-verbose.log", \
                          uflfilename[:-4]+".ffc-debug.log"]
      for filename in cleanupfilenames:
        if os.path.isfile(filename):
          os.remove(filename)

