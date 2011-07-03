import libspud
import ufltools.solver

class SpudSolver(ufltools.solver.Solver):
  def fill(self, optionpath, system):
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
