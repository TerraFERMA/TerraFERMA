import libspud
import ufltools.solver

class SpudSolver(ufltools.solver.Solver):
  """A class that stores all the information necessary to write the ufl for a system of forms (i.e. linear or bilinear) associated with a solver 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, system):
    """Fill a solver class with data describing a system of forms using libspud, the given optionpath and the system its based on."""
    self.name = libspud.get_option(optionpath+"/name")
    newoptionpath = optionpath+"/type"
    self.type = libspud.get_option(newoptionpath+"/name")
     
    if libspud.have_option(newoptionpath+"/preamble"):
      self.preamble = libspud.get_option(newoptionpath+"/preamble")+"\n"

    self.form_names = []
    self.forms = []
    self.form_symbols = []
    for i in range(libspud.option_count(newoptionpath+"/form")):
      form_optionpath = newoptionpath+"/form["+`i`+"]"
      self.form_names.append(libspud.get_option(form_optionpath+"/name"))
      self.forms.append(libspud.get_option(form_optionpath)+"\n")
      self.form_symbols.append(libspud.get_option(form_optionpath+"/ufl_symbol"))
    
    self.system = system

