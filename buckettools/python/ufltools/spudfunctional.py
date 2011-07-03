import libspud
import ufltools.functional

class SpudFunctional(ufltools.functional.Functional):
  def fill(self, optionpath, function):
    self.name     = libspud.get_option(optionpath+"/name")
    self.symbol   = libspud.get_option(optionpath+"/ufl_symbol")
    self.form     = libspud.get_option(optionpath)+"\n"
    self.function = function
  
