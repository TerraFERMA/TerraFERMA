
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
  system_ufl.append("\n")

  for j in range(libspud.option_count(system_optionpath+"/coefficient")):
    coeff_optionpath = system_optionpath+"/coefficient["+`j`+"]"

    coeff_details = get_function_details(coeff_optionpath, system_cell)

    if coeff_details["type"] == "Constant":
      system_ufl.append(usage_comment(coeff_details["name"], coeff_details["type"]))
      system_ufl.append(constant_ufl(coeff_details["symbol"], system_cell))
    else:
      system_ufl += add_element_ufl(coeff_details, system_cell)
      system_ufl.append(coefficient_ufl(coeff_details["symbol"]))
    system_ufl.append("\n")

    for k in range(libspud.option_count(coeff_optionpath+"/type/output/include_in_diagnostics/functional")):
      # Aliased coefficients don't get in here because we're not using the aliased optionpath above
      functional_optionpath = coeff_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional_details = get_diagnosticfunctional_details(functional_optionpath)
      write_diagnosticfunctional_ufl(functional_details, coeff_details, system_name, system_cell)
  
  for j in range(libspud.option_count(system_optionpath+"/nonlinear_solver")):
    solver_optionpath = system_optionpath+"/nonlinear_solver["+`j`+"]"
    solver_name       = libspud.get_option(solver_optionpath+"/name")
    solvertype_optionpath = solver_optionpath+"/type"
    solver_type           = libspud.get_option(solvertype_optionpath+"/name")

    solver_ufl = []
    solver_ufl += system_ufl
    if libspud.have_option(solvertype_optionpath+"/preamble"):
      solver_ufl.append(libspud.get_option(solvertype_optionpath+"/preamble")+"\n")

    solverform_symbols = []
    for k in range(libspud.option_count(solvertype_optionpath+"/form")):
      form_optionpath = solvertype_optionpath+"/form["+`k`+"]"
      solver_ufl.append(libspud.get_option(form_optionpath)+"\n")
      solverform_symbols.append(libspud.get_option(form_optionpath+"/ufl_symbol"))

    solver_ufl.append("\n")
    solver_ufl.append(forms_ufl(solverform_symbols))
    solver_ufl.append("\n")
    solver_ufl.append(produced_comment())

    filename = system_name+solver_name+".ufl"
    filehandle = file(filename, 'w')
    filehandle.writelines(solver_ufl)
    filehandle.close()

# and we're done
