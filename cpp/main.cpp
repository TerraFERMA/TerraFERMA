#include <dolfin.h>
#include <spud>
#include "BucketTools.h"

int main(int argc, char* argv[])
{

  buckettools::init(argc, argv);

  buckettools::SpudBucket bucket;
  bucket.fill();
  bucket.run();

  return 0;
  
}
