#include <dolfin.h>
#include <spud.h>
#include "BucketTools.h"
#include "Python.h"

int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);
  dolfin::set_log_active(true);

  const std::string options_filename = argv[1];
  Spud::load_options(options_filename);

  buckettools::SpudBucket bucket("PoissonBucket", "");
  bucket.fill();

  std::cout << bucket.str();

  return 0;
  
}
