import libspud
import ufltools.functional

class SpudFunctional(ufltools.functional.Functional):
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl) 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, function):
    """Fill a functional class with data describing that functional using libspud, the given optionpath and the function its based on."""
    self.name     = libspud.get_option(optionpath+"/name")
    self.symbol   = libspud.get_option(optionpath+"/ufl_symbol")
    self.form     = libspud.get_option(optionpath)+"\n"
    self.function = function
  