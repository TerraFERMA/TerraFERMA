
from optparse import OptionParser
import sys
import libspud
import ufltools.spudsystem 
import ufltools.spudfunction
import ufltools.spudfunctional
import ufltools.spudsolver

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file>',
                       add_help_option=True,
                       description="""This script takes a buckettools options file and writes""" +
                       """ufl files based on the nonlinear_solvers and diagnostic functionals.""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename  = argv[0]

libspud.load_options(options_filename)

for i in range(libspud.option_count("/system")):
  system = ufltools.spudsystem.SpudSystem()
  system_optionpath = "/system["+`i`+"]"
  system.fill(system_optionpath)

  for j in range(libspud.option_count(system_optionpath+"/field")):
    field_optionpath = system_optionpath+"/field["+`j`+"]"
    field = ufltools.spudfunction.SpudFunction()
    field.fill(field_optionpath, system)
     
    system.fields.append(field)

    for k in range(libspud.option_count(field_optionpath+"/type/output/include_in_diagnostics/functional")):
      functional_optionpath = field_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional = ufltools.spudfunctional.SpudFunctional()
      functional.fill(functional_optionpath, field)
      functional.write_ufl()
      del functional

  for j in range(libspud.option_count(system_optionpath+"/coefficient")):
    coeff_optionpath = system_optionpath+"/coefficient["+`j`+"]"
    coeff = ufltools.spudfunction.SpudFunction()
    coeff.fill(coeff_optionpath, system)

    system.coeffs.append(coeff)

    for k in range(libspud.option_count(coeff_optionpath+"/type/output/include_in_diagnostics/functional")):
      # Aliased coefficients don't get in here because we're not using the aliased optionpath above
      functional_optionpath = coeff_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional = ufltools.spudfunctional.SpudFunctional()
      functional.fill(functional_optionpath, coeff)
      functional.write_ufl()
      del functional
  
  for j in range(libspud.option_count(system_optionpath+"/nonlinear_solver")):
    solver_optionpath = system_optionpath+"/nonlinear_solver["+`j`+"]"
    solver = ufltools.spudsolver.SpudSolver()
    solver.fill(solver_optionpath, system)
    solver.write_ufl()
    del solver


  del system

# and we're done
