
#ifndef __USAGE_H
#define __USAGE_H


#include "Usage.h"
#include <dolfin.h>
#include <getopt.h>
#include <spud>

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // A collection of subroutines to help users at the command line
  //*****************************************************************|************************************************************//

  void usage(char *cmd);                                               // print usage to std::cerr

  void print_version(std::ostream& stream);                            // print version

  void print_environment(std::ostream& stream);                        // print relevant environment variables

  int parse_verbosity(const std::string &verbosity);                   // parse the verbosity command line argument

  void parse_arguments(int argc, char** argv);                         // parse command line arguments

}

#endif

