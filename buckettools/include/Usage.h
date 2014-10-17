// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __USAGE_H
#define __USAGE_H


#include <dolfin.h>
#include <getopt.h>
#include <spud>

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // A collection of subroutines to help users at the command line
  //*****************************************************************|************************************************************//

  void init(int argc, char** argv);                                  // parse command line arguments and initialize signal handling

  void usage(char *cmd);                                             // print usage to std::cerr

  void print_version();                                              // print version to std::cerr

  void print_environment();                                          // print relevant environment variables

  int parse_dolfin_verbosity(const std::string &verbosity);          // parse the dolfin verbosity command line argument

  int parse_verbosity(const std::string &verbosity);                 // parse the verbosity command line argument

  void parse_arguments(int argc, char** argv);                       // parse command line arguments

  std::string githash();                                             // return a string describing the git hash

  std::string compiletime();                                         // return a string describing the compilation time

}

#endif

