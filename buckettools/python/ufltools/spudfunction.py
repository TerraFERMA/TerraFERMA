import sys
import libspud
import ufltools.function

class SpudFunction(ufltools.function.Function):
  """A class that stores all the information necessary to write the ufl for a function (field or coefficient) 
     plus all the information necessary to populate that class using libspud.
     Note that the class has limited ufl production because much of this is system dependent."""

  def fill(self, optionpath, system):
    """Fill a function class with data describing that function using libspud, the given optionpath and the system its based on."""
    self.name       = libspud.get_option(optionpath+"/name")
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol")
    self.system     = system
    self.type       = libspud.get_option(optionpath+"/type/name")
    if self.type == "Aliased":
      aliasedsystem_name   = libspud.get_option(optionpath+"/type/system")
      aliasedfunction_name = libspud.get_option(optionpath+"/type/generic_function")
      
      aliasedsystem_optionpath = "/system::"+aliasedsystem_name
      if libspud.have_option(aliasedsystem_optionpath+"/field::"+aliasedfunction_name):
        aliasedfunction_optionpath = aliasedsystem_optionpath+"/field::"+aliasedfunction_name
      elif libspud.have_option(aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name):
        aliasedfunction_optionpath = aliasedsystem_optionpath+"/coefficient::"+aliasedfunction_name
      else:
        print "Unable to find aliased generic_function as either a field or a coefficient."
        sys.exit(1)
      self.type = libspud.get_option(aliasedfunction_optionpath+"/type/name")
      if self.type=="Aliased":
        print "Can't alias to an aliased function."
        sys.exit(1)
       
      # Perform a check that the meshes at least have the same name (a necessary check but is it sufficient?)
      aliasedsystemmesh_name = libspud.get_option(aliasedsystem_optionpath+"/mesh/name")
      assert(aliasedsystemmesh_name == self.system.mesh_name)
      
      newoptionpath = aliasedfunction_optionpath
    else:
      newoptionpath = optionpath

    self.rank   = libspud.get_option(newoptionpath+"/type/rank/name")
    self.family   = None
    self.degree = None
    if self.type != "Constant":
      self.family = libspud.get_option(newoptionpath+"/type/rank/element/family")
      self.degree = libspud.get_option(newoptionpath+"/type/rank/element/degree")
    
    self.size     = None
    self.shape    = None
    self.symmetry = None
    if self.rank == "Vector":
      if libspud.have_option(newoptionpath+"/type/rank/element/size"):
        self.size = libspud.get_option(newoptionpath+"/type/rank/element/size")
    elif self.rank == "Tensor":
      if libspud.have_option(newoptionpath+"/type/rank/element/shape"):
        self.shape = libspud.get_option(newoptionpath+"/type/rank/element/shape")
      if libspud.have_option(newoptionpath+"/type/rank/element/symmetry"):
        self.symmetry = True

