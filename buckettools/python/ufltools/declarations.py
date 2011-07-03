import sys
from ufltools.comments import *

def element_ufl(function, cell):
  ufl = []
  ufl.append(declaration_comment("Element", function["name"]))
  ufl_line = function["symbol"]+"_e = "
  if function["rank"] == "Scalar":
    ufl_line += "FiniteElement("
  elif function["rank"] == "Vector":
    ufl_line += "VectorElement("
  elif function["rank"] == "Tensor":
    ufl_line += "TensorElement("
  else:
    print function["rank"]
    print "Unknown rank."
    sys.exit(1)
  ufl_line += "\""+function["family"] +"\", " \
             +cell+", " \
             +`function["degree"]`
  if function["rank"] == "Vector":
    if function["size"]: ufl_line += ", size="+`function["size"]`
  elif function["rank"] == "Tensor":
    if function["shape"]: ufl_line += ", shape=("+`function["shape"][0]`+","+`function["shape"][1]`+")"
    if function["symmetry"]: ufl_line += ", symmetry=True"
  ufl_line +=")\n"
  ufl.append(ufl_line)
  ufl.append("\n")
  return ufl

def systemelement_ufl(symbol, field_symbols):
  ufl = []
  if len(field_symbols)==1:
    ufl.append(symbols_equal(symbol, field_symbols[0], suffix="_e"))
  else:
    ufl = symbol+"_e = "
    ufl += "MixedElement(["
    for s in range(len(field_symbols)-1):
      ufl += field_symbols[s]+"_e, "
    ufl += field_symbols[-1]+"_e"
    ufl += "])\n"
  return ufl

def symbols_equal(symbol_a, symbol_b, suffix=""):
  return symbol_a+suffix+" = "+symbol_b+suffix

def systemtest_ufl(symbol, field_symbols):
  ufl = []
  if len(field_symbols)==1:
    ufl.append(testfunction_ufl(symbol))
    ufl.append(testfunction_ufl(field_symbols[0]))
  else:
    ufl.append(testfunction_ufl(symbol))
    ufl.append(split_ufl(symbol, field_symbols, "_t"))
  return ufl

def systemtrial_ufl(symbol, field_symbols):
  ufl = []
  if len(field_symbols)==1:
    ufl.append(trialfunction_ufl(symbol))
    ufl.append(trialfunction_ufl(field_symbols[0]))
  else:
    ufl.append(trialfunction_ufl(symbol))
    ufl.append(split_ufl(symbol, field_symbols, "_a"))
  return ufl

def testfunction_ufl(symbol):
  return symbol+"_t = TestFunction("+symbol+"_e)\n"

def trialfunction_ufl(symbol):
  return symbol+"_a = TrialFunction("+symbol+"_e)\n"

def split_ufl(symbol, field_symbols, suffix=""):
  ufl = "("
  for s in range(len(field_symbols)-1):
    ufl += field_symbols[s]+suffix+", "
  ufl += field_symbols[-1]+suffix
  ufl += ") = split("+symbol+suffix+")\n"
  return ufl
  
def systemiterate_ufl(symbol, field_symbols):
  ufl = []
  if len(field_symbols)==1:
    ufl.append(coefficient_ufl(symbol, suffix="_i"))
    ufl.append(coefficient_ufl(field_symbols[0], suffix="_i"))
  else:
    ufl.append(coefficient_ufl(symbol, suffix="_i"))
    ufl.append(split_ufl(symbol, field_symbols, suffix="_i"))
  return ufl

def systemold_ufl(symbol, field_symbols):
  ufl = []
  if len(field_symbols)==1:
    ufl.append(coefficient_ufl(symbol, suffix="_n"))
    ufl.append(coefficient_ufl(field_symbols[0], suffix="_n"))
  else:
    ufl.append(coefficient_ufl(symbol, suffix="_n"))
    ufl.append(split_ufl(symbol, field_symbols, suffix="_n"))
  return ufl

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

def write_diagnosticfunctional_ufl(functional, function, system_name, system_cell):
  ufl = []
  if function["type"]=="Constant":
    ufl.append(declaration_comment("Constant", function["name"]))
    ufl.append(constant_ufl(function["symbol"], system_cell))
  else:
    ufl.append(declaration_comment("Element", function["name"]))
    ufl += element_ufl(function, system_cell)
    ufl.append("\n")
    ufl.append(declaration_comment("Coefficient", function["name"]))
    ufl.append(coefficient_ufl(function["symbol"]))
  ufl.append("\n")
  ufl.append(declaration_comment("Form", functional["name"]))
  ufl.append(functional["form"])
  ufl.append("\n")
  ufl.append(generic_comment("Declare non-default form names to be accessible."))
  ufl.append(forms_ufl([functional["symbol"]]))
  ufl.append("\n")
  ufl.append(produced_comment())

  filename   = system_name+function["name"]+functional["name"]+".ufl"
  filehandle = file(filename, 'w')
  filehandle.writelines(ufl)
  filehandle.close()

