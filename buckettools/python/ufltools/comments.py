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

