
from optparse import OptionParser
import sys
import libspud
from ufltools import *

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file>',
                       add_help_option=True,
                       description="""This script takes a buckettools options file and writes""" +
                       """ufl files based on the systems.""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename  = argv[0]

libspud.load_options(options_filename)

def get_function_details(optionpath, system_cell):
  function = {}
  function["name"]   = libspud.get_option(optionpath+"/name")
  function["symbol"] = libspud.get_option(optionpath+"/ufl_symbol")
  function["type"]   = libspud.get_option(optionpath+"/type/name")
  if function["type"] == "Aliased":
    aliasedsystem_name   = libspud.get_option(optionpath+"/type/system")
    aliasedfunction_name = libspud.get_option(optionpath+"/type/generic_function")
    
    aliasedsystem_optionpath = "/system::"+aliasedsystem_name
    if libspud.have_option(aliasedsystem_optionpath+"/field::"+aliasedfunction_name):
      aliasedfunction_optionpath = aliasedsystem_optionpath+"/field::"+aliasedfunction_name
    elif libspud.have_option(aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name):
      aliasedfunction_optionpath = aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name
    else:
      print "Unable to find aliased generic_function as either a field or a coefficient."
      sys.exit(1)
    function["type"] = libspud.get_option(aliasedfunction_optionpath+"/type/name")
    if function["type"]=="Aliased":
      print "Can't alias to an aliased function."
      sys.exit(1)
     
    # Perform a check that the meshes at least have the same cell (a necessary but not sufficient check)
    aliasedsystemmesh_name       = libspud.get_option(aliasedsystem_optionpath+"/mesh/name")
    aliasedsystemmesh_optionpath = "/geometry/mesh::"+aliasedsystemmesh_name
    aliasedsystem_cell           = libspud.get_option(aliasedsystemmesh_optionpath+"/cell")
    assert(aliasedsystem_cell == system_cell)
    
    newoptionpath = aliasedfunction_optionpath
  else:
    newoptionpath = optionpath

  function["optionpath"] = optionpath
  function["rank"]   = libspud.get_option(newoptionpath+"/type/rank/name")
  function["family"]   = None
  function["degree"] = None
  if function["type"] != "Constant":
    function["family"] = libspud.get_option(newoptionpath+"/type/rank/element/family")
    function["degree"] = libspud.get_option(newoptionpath+"/type/rank/element/degree")
  
  function["size"]     = None
  function["shape"]    = None
  function["symmetry"] = None
  if function["rank"] == "Vector":
    if libspud.have_option(newoptionpath+"/type/rank/element/size"):
      function["size"] = libspud.get_option(newoptionpath+"/type/rank/element/size")
  elif function["rank"] == "Tensor":
    if libspud.have_option(newoptionpath+"/type/rank/element/shape"):
      function["shape"] = libspud.get_option(newoptionpath+"/type/rank/element/shape")
    if libspud.have_option(newoptionpath+"/type/rank/element/symmetry"):
      function["symmetry"] = True

  return function

def get_diagnosticfunctional_details(optionpath):
  functional = {}
  functional["name"]   = libspud.get_option(optionpath+"/name")
  functional["symbol"] = libspud.get_option(optionpath+"/ufl_symbol")
  functional["form"]   = libspud.get_option(optionpath)+"\n"
  return functional
  
for i in range(libspud.option_count("/system")):
  system_optionpath = "/system["+`i`+"]"
  system_name       = libspud.get_option(system_optionpath+"/name")
  system_symbol     = libspud.get_option(system_optionpath+"/ufl_symbol")

  systemmesh_name       = libspud.get_option(system_optionpath+"/mesh/name")
  systemmesh_optionpath = "/geometry/mesh::"+systemmesh_name
  system_cell           = libspud.get_option(systemmesh_optionpath+"/cell")

  systemfield_symbols = []
  systemcoeff_symbols = []
  systemconst_symbols = []
  system_ufl          = []

  for j in range(libspud.option_count(system_optionpath+"/field")):
    field_optionpath = system_optionpath+"/field["+`j`+"]"
     
    field_details = get_function_details(field_optionpath, system_cell)    
    
    system_ufl += add_element_ufl(field_details, system_cell)
    
    systemfield_symbols.append(field_details["symbol"])

    for k in range(libspud.option_count(field_optionpath+"/type/output/include_in_diagnostics/functional")):
      functional_optionpath = field_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional_details = get_diagnosticfunctional_details(functional_optionpath)
      write_diagnosticfunctional_ufl(functional_details, field_details, system_name, system_cell)

  for j in range(libspud.option_count(system_optionpath+"/coefficient")):
    coeff_optionpath = system_optionpath+"/coefficient["+`j`+"]"

    coeff_details = get_function_details(coeff_optionpath, system_cell)

    system_ufl += add_element_ufl(coeff_details, system_cell)

    if coeff_details["type"]=="Constant":
      systemconst_symbols.append(coeff_details["symbol"])
    else:
      systemcoeff_symbols.append(coeff_details["symbol"])
    
    for k in range(libspud.option_count(coeff_optionpath+"/type/output/include_in_diagnostics/functional")):
      # Aliased coefficients don't get in here because we're not using the aliased optionpath above
      functional_optionpath = coeff_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional_details = get_diagnosticfunctional_details(functional_optionpath)
      write_diagnosticfunctional_ufl(functional_details, coeff_details, system_name, system_cell)
  

  system_ufl.append("\n")
  system_ufl.append(declaration_comment("Element", system_name))
  system_ufl.append(systemelement_ufl(system_symbol, systemfield_symbols))
  system_ufl.append("\n")
  system_ufl += systemtest_ufl(system_symbol, systemfield_symbols)
  system_ufl.append("\n")
  system_ufl += systemtrial_ufl(system_symbol, systemfield_symbols)
  system_ufl.append("\n")
  system_ufl += systemiterate_ufl(system_symbol, systemfield_symbols)
  system_ufl.append("\n")
  system_ufl += systemold_ufl(system_symbol, systemfield_symbols)

  filename = system_name+".ufl"
  filehandle = file(filename, 'w')
  filehandle.writelines(system_ufl)
  filehandle.close()

  #for j in range(libspud.option_count("/system["+`i`+"]/nonlinear_solver")):
  #  solver_optionpath = system_optionpath+"/nonlinear_solver["+`j`+"]"
  #  solver_name       = libspud.get_option(solver_optionpath+"/name")
  #  solver_type       = libspud.get_option(solver_optionpath+"/type/name")
  #  solver_optionpath += 
  
# and we're done
