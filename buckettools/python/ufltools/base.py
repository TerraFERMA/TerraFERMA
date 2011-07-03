import sys

# Functions used to provide comments in ufls
def comment(string):
  return "# "+string+"\n"

def declaration_comment(thing, type, name):
  return comment(thing+" declaration for "+type+": "+name)

def produced_comment():
  return comment("Produced by: "+" ".join(sys.argv))

def equal_ufl(symbol_a, symbol_b, suffix=""):
  return symbol_a+suffix+" = "+symbol_b+suffix

def coefficient_ufl(symbol, suffix=""):
  return symbol+suffix+" = Coefficient("+symbol+"_e)\n"

def constant_ufl(symbol, cell, suffix=""):
  return symbol+suffix+" = Constant("+cell+")\n"

def forms_ufl(form_symbols):
  ufl = "forms = ["
  for s in range(len(form_symbols)-1):
    ufl += form_symbols[s]+", "
  ufl += form_symbols[-1]
  ufl += "]\n"
  return ufl

def testfunction_ufl(symbol):
  return symbol+"_t = TestFunction("+symbol+"_e)\n"

def trialfunction_ufl(symbol):
  return symbol+"_a = TrialFunction("+symbol+"_e)\n"

