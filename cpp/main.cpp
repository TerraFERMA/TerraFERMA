#include <dolfin.h>
#include <spud>
#include "BucketTools.h"
#include "Python.h"

int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);
  dolfin::set_log_active(true);
  dolfin::set_log_level(dolfin::DBG);

  std::stringstream buffer;
  Spud::OptionError serr;

  const std::string options_filename = argv[1];
  Spud::load_options(options_filename);

  std::string basename;
  buffer.str(""); buffer << "/io/output_base_name";
  serr = Spud::get_option(buffer.str(), basename); buckettools::spud_err(buffer.str(), serr);

  buckettools::SpudBucket bucket(basename, "");

  bucket.fill();
  bucket.output();

  buckettools::DiagnosticsFile diagfile(basename+".stat");
  diagfile.write_header(bucket, false);
  diagfile.write_data(bucket);

  bucket.run();

  bucket.output();

  diagfile.write_data(bucket);

  std::cout << bucket.str();

  return 0;
  
}
