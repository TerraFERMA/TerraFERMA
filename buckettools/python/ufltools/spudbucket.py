import libspud
import ufltools.bucket
import ufltools.spud

class SpudBucket(ufltools.bucket.Bucket):
  def fill(self):
    """Fill a bucket class with data describing a set of mixedfunctionspace systems using libspud, the given optionpath."""

    self.systems = []

    # loop over the systems in the options tree
    for i in range(libspud.option_count("/system")):
      system_optionpath = "/system["+`i`+"]"
      system = ufltools.spud.SpudSystem()
      # get all the information about this system from the options dictionary
      system.fill(system_optionpath, self)
      # let the bucket know about this system
      self.systems.append(system)
      # done with this system
      del system

