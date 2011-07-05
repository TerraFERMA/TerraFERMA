
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

cpp       = []
cpp.append("#ifndef __SYSTEMS_WRAPPER_H\n")
cpp.append("#define __SYSTEMS_WRAPPER_H\n")
cpp.append("\n")
cpp.append("#include \"Bucket.h\"\n")
cpp.append("#include <dolfin.h>\n")
cpp.append("\n")

include_cpp = []

functionspace_cpp         = []
functionspace_cpp.append("  // A function to return a functionspace from a system given a mesh and a solvername.\n")
functionspace_cpp.append("  FunctionSpace_ptr fetch_functionspace(std::string systemname, std::string solvername, Mesh_ptr mesh)\n")
functionspace_cpp.append("  {\n")
functionspace_cpp.append("    switch(systemname)\n")
functionspace_cpp.append("    {\n")

coefficientspace_cpp         = []
coefficientspace_cpp.append("  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a functionname.\n")
coefficientspace_cpp.append("  FunctionSpace_ptr fetch_coefficientspace(std::string systemname, std::string solvername, std::string functionname, Mesh_ptr mesh)\n")
coefficientspace_cpp.append("  {\n")
coefficientspace_cpp.append("    switch(systemname)\n")
coefficientspace_cpp.append("    {\n")

functional_cpp            = []
functional_cpp.append("  // A function to return a functional from a system-function set given a mesh and a functionalname.\n")
functional_cpp.append("  Form_ptr fetch_functional(std::string systemname, std::string functionname, std::string functionalname, Mesh_ptr mesh)\n")
functional_cpp.append("  {\n")
functional_cpp.append("    switch(systemname)\n")
functional_cpp.append("    {\n")

form_cpp = []
form_cpp.append("  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.\n")
form_cpp.append("  Form_ptr fetch_form(std::string systemname, std::string solvername, std::string solvertype, std::string formname, FunctionSpace_ptr functionspace)\n")
form_cpp.append("  {\n")
form_cpp.append("    switch(systemname)\n")
form_cpp.append("    {\n")

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
      # let the field know about this functional
      field.functionals.append(functional)
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
      # let the coefficient know about this functional
      coeff.functionals.append(functional)
      # done with this functional
      del functional
    # remove the local copy of this coefficient
    del coeff
  for j in range(libspud.option_count(system_optionpath+"/nonlinear_solver")):
    solver_optionpath = system_optionpath+"/nonlinear_solver["+`j`+"]"
    solver = ufltools.spud.SpudSolver()
    # get all the information about this nonlinear solver from the options dictionary
    solver.fill(solver_optionpath, system)
    # let the system know about this solver
    system.solvers.append(solver)
    # done with this nonlinear solver
    del solver
  include_cpp    += system.include_cpp()
  functionspace_cpp += system.functionspace_cpp()
  coefficientspace_cpp += system.coefficientspace_cpp()
  functional_cpp += system.functional_cpp()
  form_cpp += system.form_cpp()
  # done with this system
  del system

functionspace_cpp.append("      default:\n")
functionspace_cpp.append("        dolfin::error(\"Unknown systemname in fetch_functionspace\");\n")
functionspace_cpp.append("    }\n")
functionspace_cpp.append("    return functionspace;\n")
functionspace_cpp.append("  }\n")

coefficientspace_cpp.append("      default:\n")
coefficientspace_cpp.append("        dolfin::error(\"Unknown systemname in fetch_coefficientspace\");\n")
coefficientspace_cpp.append("    }\n")
coefficientspace_cpp.append("    return coefficientspace;\n")
coefficientspace_cpp.append("  }\n")

functional_cpp.append("      default:\n")
functional_cpp.append("        dolfin::error(\"Unknown systemname in fetch_functional\");\n")
functional_cpp.append("    }\n")
functional_cpp.append("    return functional;\n")
functional_cpp.append("  }\n")

form_cpp.append("      default:\n")
form_cpp.append("        dolfin::error(\"Unknown systemname in fetch_form\");\n")
form_cpp.append("    }\n")
form_cpp.append("    return form;\n")
form_cpp.append("  }\n")

cpp += include_cpp
cpp.append("\n")
cpp.append("namespace buckettools\n")
cpp.append("{\n")
cpp += functionspace_cpp
cpp.append("\n")
cpp += coefficientspace_cpp
cpp.append("\n")
cpp += form_cpp
cpp.append("\n")
cpp += functional_cpp
cpp.append("}\n")
cpp.append("\n")
cpp.append("#endif\n")

filehandle = file("SystemsWrapper.h", 'w')
filehandle.writelines(cpp)
filehandle.close()

# and we're done

