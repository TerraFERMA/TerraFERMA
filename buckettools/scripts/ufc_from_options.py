
from optparse import OptionParser
import sys
import libspud
import ufltools.spud

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file>',
                       add_help_option=True,
                       description="""This script takes a buckettools options file and writes """ +
                       """ufl files based on the nonlinear_solver and diagnostic functional options.  """ +
                       """It also provides a cpp header file wrapping the namespaces of the ufc """ +
                       """corresponding to these ufl files.""" )

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename  = argv[0]

# load the options tree
libspud.load_options(options_filename)

bucket = ufltools.spud.SpudBucket()
# populate the bucket based on the options file loaded above
bucket.fill()
# write out the ufl files described by the options tree
bucket.write_ufc()
# write a cpp header file to wrap the namespaces of the corresponding ufc
bucket.write_systemfunctionals_cpp()
bucket.write_systemsolvers_cpp()

# and we're done

