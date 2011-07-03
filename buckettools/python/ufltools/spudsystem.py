import ufltools.system
import libspud

class SpudSystem(ufltools.system.System):
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space) 
     plus all the information necessary to populate that class using libspud."""
  
  def fill(self, optionpath):
    """Fill a system class with data describing that system using libspud and the given optionpath."""
    self.name       = libspud.get_option(optionpath+"/name")
    self.optionpath = optionpath
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol")

    self.mesh_name        = libspud.get_option(optionpath+"/mesh/name")
    mesh_optionpath  = "/geometry/mesh::"+self.mesh_name
    self.cell             = libspud.get_option(mesh_optionpath+"/cell")

    self.fields = []
    self.coeffs = []
    
 
