
from optparse import OptionParser
import sys
import libspud
import ufltools.spud

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file>',
                       add_help_option=True,
                       description="""This script takes a buckettools options file and writes""" +
                       """ufl files based on the nonlinear_solver and diagnostic functional options.""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename  = argv[0]

# load the options tree
libspud.load_options(options_filename)

# loop over the systems in the options tree
for i in range(libspud.option_count("/system")):
  system_optionpath = "/system["+`i`+"]"
  system = ufltools.spud.SpudSystem()
  # get all the information about this system from the options dictionary
  system.fill(system_optionpath)

  for j in range(libspud.option_count(system_optionpath+"/field")):
    field_optionpath = system_optionpath+"/field["+`j`+"]"
    field = ufltools.spud.SpudFunction()
    # get all the information about this field from the options dictionary
    field.fill(field_optionpath, system)
    # let the system know about this field
    system.fields.append(field)

    for k in range(libspud.option_count(field_optionpath+"/type/output/include_in_diagnostics/functional")):
      functional_optionpath = field_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional = ufltools.spud.SpudFunctional()
      # get all the information about this functional from the options dictionary
      functional.fill(functional_optionpath, field)
      # write some ufl to disk
      functional.write_ufl()
      # done with this functional
      del functional

    # remove the local copy of this field
    del field

  for j in range(libspud.option_count(system_optionpath+"/coefficient")):
    coeff_optionpath = system_optionpath+"/coefficient["+`j`+"]"
    coeff = ufltools.spud.SpudFunction()
    # get all the information about this coefficient from the options dictionary
    coeff.fill(coeff_optionpath, system)
    # let the system know about this coefficient
    system.coeffs.append(coeff)

    for k in range(libspud.option_count(coeff_optionpath+"/type/output/include_in_diagnostics/functional")):
      # aliased coefficients don't get in here because we're not using the aliased optionpath in the line above
      functional_optionpath = coeff_optionpath+"/type/output/include_in_diagnostics/functional["+`k`+"]"
      functional = ufltools.spud.SpudFunctional()
      # get all the information about this functional from the options dictionary
      functional.fill(functional_optionpath, coeff)
      # write some ufl to disk
      functional.write_ufl()
      # done with this functional
      del functional

    # remove the local copy of this coefficient
    del coeff
  
  for j in range(libspud.option_count(system_optionpath+"/nonlinear_solver")):
    solver_optionpath = system_optionpath+"/nonlinear_solver["+`j`+"]"
    solver = ufltools.spud.SpudSolver()
    # get all the information about this nonlinear solver from the options dictionary
    solver.fill(solver_optionpath, system)
    # write some ufl to disk
    solver.write_ufl()
    # done with this nonlinear solver
    del solver

  # done with this system
  del system

# and we're done
