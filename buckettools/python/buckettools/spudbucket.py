import libspud
import buckettools.bucket
import buckettools.spud

class SpudBucket(buckettools.bucket.Bucket):
  def fill(self):
    """Fill a bucket class with data describing a set of mixedfunctionspace systems using libspud, the given optionpath."""

    self.meshes = {}
    # loop over the meshes in the options tree
    for i in range(libspud.option_count("/geometry/mesh")):
      mesh_optionpath = "/geometry/mesh["+`i`+"]"
      mesh_name = libspud.get_option(mesh_optionpath+"/name")
      self.meshes[mesh_name] = libspud.get_option(mesh_optionpath+"/source/cell")

    visualization_optionpath = "/io/visualization/element"
    self.viselementfamily = libspud.get_option(visualization_optionpath+"/family")
    self.viselementdegree = libspud.get_option(visualization_optionpath+"/degree")

    parameters_optionpath = "/global_parameters/ufl"
    if libspud.have_option(parameters_optionpath):
      self.parameters = libspud.get_option(parameters_optionpath)

    self.systems = []
    # loop over the systems in the options tree
    for i in range(libspud.option_count("/system")):
      system_optionpath = "/system["+`i`+"]"
      system = buckettools.spud.SpudSystemBucket()
      # get all the information about this system from the options dictionary
      system.fill(system_optionpath, self)
      # let the bucket know about this system
      self.systems.append(system)
      # done with this system
      del system

