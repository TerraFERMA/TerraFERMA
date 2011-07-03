import sys

# Functions used to provide comments in ufls
def generic_comment(comment):
  return "# "+comment+"\n"

def declaration_comment(type, name):
  return generic_comment(type+" declaration for "+name)

def usage_comment(name, usage):
  lvowels = ['a', 'i', 'e', 'o', 'u']
  uvowels = [l.upper() for l in lvowels] 
  if usage[0] in lvowels or usage[0] in uvowels:
    return generic_comment(name+" will be used as an "+usage)
  else:
    return generic_comment(name+" will be used as a "+usage)

def produced_comment():
  return generic_comment("Produced by: "+" ".join(sys.argv))

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

