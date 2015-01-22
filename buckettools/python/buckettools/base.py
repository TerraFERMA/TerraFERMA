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

# Functions used to provide comments in ufls
def comment(string):
  """Returns a ufl formatted comment based around the provided string."""
  return "# "+string+"\n"

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
  return symbol+suffix+" = Coefficient("+symbol+"_e)\n"

def forms_ufl(form_symbols):
  """Returns a ufl string declaring a set of forms with potentially non-default symbols."""
  ufl = "forms = ["
  for s in range(len(form_symbols)-1):
    ufl += form_symbols[s]+", "
  ufl += form_symbols[-1]
  ufl += "]\n"
  return ufl

def testfunction_ufl(symbol):
  """Returns a ufl string declaring a test function on an element."""
  return symbol+"_t = TestFunction("+symbol+"_e)\n"

def trialfunction_ufl(symbol):
  """Returns a ufl string declaring a trial function on an element."""
  return symbol+"_a = TrialFunction("+symbol+"_e)\n"

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

def form_symbols(form):
  return [symbol for line in form.split("\n") if not line.lstrip().startswith("#") for symbol in re.findall(r"\b\w+?\b", line, re.I) if symbol not in ufl_reserved()]

def forms_symbols(forms):
  symbols = []
  for form in forms:
    symbols += form_symbols(form)
  return symbols

def ffc(namespace, quadrature_rule, quadrature_degree):
  uflfilename = namespace+".ufl"

  try:
    checksum = hashlib.md5(open(uflfilename).read()).hexdigest()
  except:
    checksum = None

  try:
    headerfile = open(namespace+".h")
  except IOError:
    headerfile = None

  if headerfile:
    rebuild = checksum != hashlib.md5(open(uflfilename+".temp").read()).hexdigest()

    headertext = headerfile.read()
    
    qdi  = headertext.find("quadrature_degree")
    qdin = headertext.find("\n", qdi, -1)
    if (qdi != -1) and (qdin != -1):
      qd_old = headertext[qdi:qdin].split(" ")[-1]
      qd_new = "'"+`quadrature_degree`+"'" if quadrature_degree is not None else "'auto'"
      rebuild = rebuild or qd_new != qd_old
    else: # if we don't know what the old quadrature degree was we always want to (re)build
      rebuild = True

    qri  = headertext.find("quadrature_rule")
    qrin = headertext.find("\n", qri, -1)
    if (qri != -1) and (qrin != -1):
      qr_old = headertext[qri:qrin].split(" ")[-1]
      qr_new = "'"+quadrature_rule+"'" 
      rebuild = rebuild or qr_new != qr_old
    else: # if we don't know what the old quadrature rule was we always want to (re)build
      rebuild = True
    
  else:   # if we don't have a header file we always want to (re)build
    rebuild = True

  if rebuild:
    # files and/or quadrature degree have changed
    shutil.copy(uflfilename+".temp", uflfilename)
    command = ["ffc", "-l", "dolfin", "-O", "-r", "quadrature"]
    if quadrature_degree is not None:
      command += ["-fquadrature_degree="+`quadrature_degree`]
    command += ["-fsplit"]
    command += ["-q", quadrature_rule]
    command += [uflfilename]
    try:
      subprocess.check_call(command)
    except:
      print "ERROR while calling ffc on file ", uflfilename
      sys.exit(1)

