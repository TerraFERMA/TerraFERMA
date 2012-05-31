import buckettools.systembucket
import buckettools.spud
import libspud

class SpudSystemBucket(buckettools.systembucket.SystemBucket):
  """A class that stores all the information necessary to write the ufl for a system (i.e. mixed function space) 
     plus all the information necessary to populate that class using libspud."""
  
  def fill(self, optionpath, bucket):
    """Fill a system class with data describing that system using libspud and the given optionpath."""
    self.name       = libspud.get_option(optionpath+"/name")
    self.optionpath = optionpath
    self.symbol     = libspud.get_option(optionpath+"/ufl_symbol").split("\n")[0]
    self.bucket     = bucket

    self.mesh_name        = libspud.get_option(optionpath+"/mesh/name")
    mesh_optionpath  = "/geometry/mesh::"+self.mesh_name
    self.cell             = libspud.get_option(mesh_optionpath+"/source/cell")

    self.fields = []
    for j in range(libspud.option_count(optionpath+"/field")):
      field_optionpath = optionpath+"/field["+`j`+"]"
      field = buckettools.spud.SpudFunctionBucket()
      # get all the information about this field from the options dictionary
      field.fill(field_optionpath, self, j)
      # let the system know about this field
      self.fields.append(field)
      # remove the local copy of this field
      del field

    self.coeffs = []
    for j in range(libspud.option_count(optionpath+"/coefficient")):
      coeff_optionpath = optionpath+"/coefficient["+`j`+"]"
      coeff = buckettools.spud.SpudFunctionBucket()
      # get all the information about this coefficient from the options dictionary
      coeff.fill(coeff_optionpath, self, j)
      # let the system know about this coefficient
      self.coeffs.append(coeff)
      # remove the local copy of this coefficient
      del coeff
    
    self.special_coeffs = []
    if libspud.have_option("/timestepping"):
      coeff_optionpath = "/timestepping/timestep/coefficient::Timestep"
      coeff = buckettools.spud.SpudFunctionBucket()
      # get all the information about this coefficient from the options dictionary
      coeff.fill(coeff_optionpath, self, 0)
      # let the system know about this coefficient
      self.special_coeffs.append(coeff)
      # remove the local copy of this coefficient
      del coeff
    
    self.solvers = []
    for j in range(libspud.option_count(optionpath+"/nonlinear_solver")):
      solver_optionpath = optionpath+"/nonlinear_solver["+`j`+"]"
      solver = buckettools.spud.SpudSolverBucket()
      # get all the information about this nonlinear solver from the options dictionary
      solver.fill(solver_optionpath, self)
      # let the system know about this solver
      self.solvers.append(solver)
      # done with this nonlinear solver
      del solver

