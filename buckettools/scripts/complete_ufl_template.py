
from optparse import OptionParser
import sys
from lxml import etree
from string import Template as template

# Let's start by parsing the options
optparser=OptionParser(usage='usage: %prog <options-file> <template-file> <ufl-file-name>',
                       add_help_option=True,
                       description="""This script takes an options file and a template ufl """ +
                       """and converts the template into a full ufl file based on a """ +
                       """dictionary constructed from the placeholder attributes in """ +
                       """the options file.""")

(options, argv) = optparser.parse_args()

if len(argv)<3:
    optparser.print_help()
    sys.exit(1)

# the options file name
options_filename = argv[0]
# the template file name
template_filename = argv[1]
# the output ufl file name
ufl_filename = argv[2]

# load the options using element tree
tree = etree.parse(options_filename)

# these element tags are reserved by diamond (and spud)
spud_reserved = ['comment', 'string_value', 'integer_value', 'real_value']

# compile a dictionary mapping placeholder
# to the value in the tree
listing = []
placeholders = tree.findall('.//*[@placeholder]')
for placeholder in placeholders:
  index = placeholder.attrib['placeholder']
  # get the value from the placeholder
  value = placeholder.text
  # but if we have a spud reserved name underneath
  # then replace the value from that
  for child in placeholder.getchildren():
    if child.tag in spud_reserved:
      value = child.text
      break
  if value == None:
    print 'ERROR: no value for placeholder id:', index
    sys.exit(1)
  listing.append((index, value))

dictionary = dict(listing)

# open the template and output ufl files
tempfile = file(template_filename, 'r')
uflfile = file(ufl_filename, 'w')

# use python templating to convert the
# placeholders indices to their values
ufllines = []
templines = tempfile.readlines()
for templine in templines:
  uflline = template(templine).substitute(dictionary)
  ufllines.append(uflline)

# write the ufl file
uflfile.writelines(ufllines)
# and we're done
