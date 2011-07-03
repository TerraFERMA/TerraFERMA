
from optparse import OptionParser
import sys
import libspud
from string import Template as template

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file> <template-header-file> <header-file>',
                       add_help_option=True,
                       description="""This script takes an options file and a template header file""" +
                       """and converts the template into a full header.""")

(options, argv) = optparser.parse_args()

if len(argv)<3:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename  = argv[0]
# the template header file name
template_filename = argv[1]
# the output header file name
header_filename   = argv[2]

include_ufcs       = ""
case_functionspace = ""
case_form_0        = ""
case_form_1        = ""
case_form_2        = ""


libspud.load_options(options_filename)

for i in range(libspud.option_count("/system")):
  namespace = libspud.get_option("/system["+`i`+"]/name")
  
  include_ufcs += "#include \""+namespace+".h\"\n"
  
  case_functionspace += "      case \""+namespace+"\":\n"
  case_functionspace += "        FunctionSpace_ptr sysspace(new "+namespace+"::FunctionSpace(*mesh));\n"
  case_functionspace += "        break;\n\n"
  
  case_form_0 += "      case \""+namespace+"\":\n"
  case_form_0 += "        Form_ptr form_0(new "+namespace+"::Form_0(*fs, *fs));\n"
  case_form_0 += "        break;\n\n"

  case_form_1 += "      case \""+namespace+"\":\n"
  case_form_1 += "        Form_ptr form_1(new "+namespace+"::Form_0(*fs));\n"
  case_form_1 += "        break;\n\n"

  case_form_2 += "      case \""+namespace+"\":\n"
  case_form_2 += "        Form_ptr form_2(new "+namespace+"::Form_0(*fs));\n"
  case_form_2 += "        break;\n\n"

dictionary = {                                           \
               'include_ufcs'      : include_ufcs,       \
               'case_functionspace': case_functionspace, \
               'case_form_0'       : case_form_0,        \
               'case_form_1'       : case_form_1,        \
               'case_form_2'       : case_form_2,        \
             }

# open the template and output ufl files
tempfile = file(template_filename, 'r')
templines = tempfile.read()
tempfile.close()

headlines = template(templines).substitute(dictionary)

# write the header file
headfile = file(header_filename, 'w')
headfile.writelines(headlines)
headfile.close()
# and we're done
