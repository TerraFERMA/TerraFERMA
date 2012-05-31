import libspud
import ufltools.functionalbucket

class SpudFunctionalBucket(ufltools.functionalbucket.FunctionalBucket):
  """A class that stores all the information necessary to write the ufl for a functional (i.e. scalar valued returning ufl) 
     plus all the information necessary to populate that class using libspud."""

  def fill(self, optionpath, function):
    """Fill a functional class with data describing that functional using libspud, the given optionpath and the function its based on."""
    try:
      self.name   = libspud.get_option(optionpath+"/name")
    except libspud.SpudKeyError:
      self.name   = ""
    self.symbol   = libspud.get_option(optionpath+"/ufl_symbol").split("\n")[0]
    self.form     = libspud.get_option(optionpath)+"\n"
    self.function = function
    if libspud.have_option(optionpath+"/quadrature_degree"):
      self.quadrature_degree = libspud.get_option(optionpath+"/quadrature_degree")


  
