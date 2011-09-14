#include <dolfin.h>
#include <spud>
#include "BucketTools.h"

int main(int argc, char* argv[])
{
  buckettools::parse_arguments(argc,argv);

  buckettools::SpudBucket bucket;

  bucket.fill();

  bucket.run();

  std::cout << bucket.str();

  return 0;
  
}
