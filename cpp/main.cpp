#include <dolfin.h>
#include <spud>
#include "BucketTools.h"
#include "Python.h"

int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);
  dolfin::set_log_active(true);
  dolfin::set_log_level(dolfin::DBG);

  const std::string options_filename = argv[1];
  Spud::load_options(options_filename);

  buckettools::SpudBucket bucket;

  bucket.fill();

  bucket.run();

  std::cout << bucket.str();

  return 0;
  
}
