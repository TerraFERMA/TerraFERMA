import sys

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

def constant_ufl(symbol, cell, suffix=""):
  """Returns a ufl string declaring a constant coefficient on a cell."""
  return symbol+suffix+" = Constant("+cell+")\n"

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

