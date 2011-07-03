import ufltools.system
import libspud

class SpudSystem(ufltools.system.System):
  def fill(self, optionpath):
    self.name       = libspud.get_option(optionpath+"/name")
    self.optionpath = optionpath
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol")

    self.mesh_name        = libspud.get_option(optionpath+"/mesh/name")
    mesh_optionpath  = "/geometry/mesh::"+self.mesh_name
    self.cell             = libspud.get_option(mesh_optionpath+"/cell")

    self.fields = []
    self.coeffs = []
    
 
