import sys
import re

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
  return ["", "_i", "_n"]

def form_symbols(form):
  return [symbol for line in form.split("\n") if not line.lstrip().startswith("#") for symbol in re.findall(r"\b\w+?\b", line, re.I) if symbol not in ufl_reserved()]

def forms_symbols(forms):
  symbols = []
  for form in forms:
    symbols += form_symbols(form)
  return symbols


