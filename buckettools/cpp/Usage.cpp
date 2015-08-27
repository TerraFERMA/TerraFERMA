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


#include "Usage.h"
#include "builddefs.h"
#include <dolfin.h>
#include <getopt.h>
#include <spud>
#include "SignalHandler.h"
#include "EventHandler.h"
#include "SigIntEventHandler.h"
#include <signal.h>
#include <cstdio>
#include "Logger.h"
#include "BucketPETScBase.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// parse the command line arguments and set up the signal handler
//*******************************************************************|************************************************************//
void buckettools::init(int argc, char** argv)
{

  parse_arguments(argc,argv);

  SigIntEventHandler_ptr sigint_handler( new SigIntEventHandler );
  (*SignalHandler::instance()).register_handler(SIGINT, sigint_handler);

}

//*******************************************************************|************************************************************//
// print the recommended usage
//*******************************************************************|************************************************************//
void buckettools::usage(char *cmd)
{
  print_version();
  std::stringstream s; s.str("");
  s << std::endl << "Usage: " << cmd << " [options ...] <simulation-file>" << std::endl
      << std::endl << "Options:" << std::endl
      <<" -v <level>, --verbose <level>" << std::endl << "\tVerbose output, defaults to WARNING (30) if unset. " << std::endl
      <<"\tAvailable options: ERROR (40), WARNING (30), INFO (20), DEBUG (10), DBG (10) or any integer." << std::endl
      <<" -d <level>, --dolfin-verbose <level>" << std::endl << "\tVerbose DOLFIN output, defaults to match -v if unset." << std::endl
      <<"\tAvailable options: CRITICAL (50), ERROR (40), WARNING (30), INFO (20), PROGRESS (16), TRACE (13), DEBUG (10), DBG (10) or any integer." << std::endl
      <<" -p, --petsc-info" << std::endl << "\tVerbose PETSc output." << std::endl
      <<" -l, --log" << std::endl << "\tCreate log (redirects stdout) and error (redirects stderr) files for each process." << std::endl
      <<" -V, --version" << std::endl << "\tPrints version information then exits." << std::endl
      <<" -h, --help" << std::endl << "\tHelp! Prints this message then exits.";
  log(ERROR, s.str());
  return;
}

//*******************************************************************|************************************************************//
// print the version number
//*******************************************************************|************************************************************//
void buckettools::print_version()
{
  std::stringstream s; s.str("");
  s << "GitHash: " << githash() << std::endl
      <<"CompileTime: "<< compiletime();
  log(ERROR, s.str());
  return;
}

//*******************************************************************|************************************************************//
// print relevant environment variables to an std::ostream
//*******************************************************************|************************************************************//
void buckettools::print_environment()
{

  const char *relevant_variables[]={"PETSC_OPTIONS"};

  int no_relevant_variables=1;
  int variables_found=0;

  std::stringstream s; s.str("");
  s << "Relevant environment variables:" << std::endl;
  for(int i=0; i<no_relevant_variables; i++) {
    char *env=getenv(relevant_variables[i]);
    if(env!=NULL) {
      s << "  " << relevant_variables[i]<<" = " << env <<std::endl;
      variables_found++;
    }
  };
  if (!variables_found)
  {
    s << "  none" << std::endl;
  }
  log(INFO, s.str());
  return;
}

//*******************************************************************|************************************************************//
// parse the command line argument for verbosity to an integer
//*******************************************************************|************************************************************//
int buckettools::parse_verbosity(const std::string &verbosity)
{
  int v;
  if (verbosity=="ERROR")
  {
    v = ERROR;
  }
  else if (verbosity=="WARNING")
  {
    v = WARNING;
  }
  else if (verbosity=="INFO")
  {
    v = INFO;
  }
  else if (verbosity=="DEBUG")
  {
    v = DBG;
  }
  else if (verbosity=="DBG")
  {
    v = DBG;
  }
  else
  {
    v = atoi(verbosity.c_str());
  }
  return v;
}

//*******************************************************************|************************************************************//
// parse the command line argument for dolfin verbosity to an integer
//*******************************************************************|************************************************************//
int buckettools::parse_dolfin_verbosity(const std::string &verbosity)
{
  int v;
  if (verbosity=="CRITICAL")
  {
    v = dolfin::CRITICAL;
  }
  else if (verbosity=="ERROR")
  {
    v = dolfin::ERROR;
  }
  else if (verbosity=="WARNING")
  {
    v = dolfin::WARNING;
  }
  else if (verbosity=="INFO")
  {
    v = dolfin::INFO;
  }
  else if (verbosity=="PROGRESS")
  {
    v = dolfin::PROGRESS;
  }
  else if (verbosity=="TRACE")
  {
    v = dolfin::TRACE;
  }
  else if (verbosity=="DEBUG")
  {
    v = dolfin::DBG;
  }
  else if (verbosity=="DBG")
  {
    v = dolfin::DBG;
  }
  else
  {
    v = atoi(verbosity.c_str());
  }
  return v;
}

//*******************************************************************|************************************************************//
// parse the command line arguments and set various options based on them
//*******************************************************************|************************************************************//
void buckettools::parse_arguments(int argc, char** argv)
{
  struct option long_options[] = {                                   // a structure linking long option names with their short equivalents
    {"help",           no_argument,       0, 'h'},
    {"log",            no_argument,       0, 'l'},
    {"petsc-info",     no_argument,       0, 'p'},
    {"verbose",        required_argument, 0, 'v'},
    {"dolfin-verbose", required_argument, 0, 'd'},
    {"version",        no_argument,       0, 'V'},
    {0,                0,                 0, 0}                      // terminated with an array of zeros
  };
  int option_index = 0;
  int verbosity, dolfinverbosity;
  int c;
  std::map< std::string, std::string > command_line_options;

  int petscargc;                                                    // dummy argc and argv to stop our command line arguments being
  char** petscargv;                                                 // interpretted by petsc usage.  FIXME: ugly hack!
  petscargc = 1;
  petscargv = new char*[1];
  petscargv[0] = new char[strlen(argv[0]) + 1];
  strncpy(petscargv[0], argv[0], strlen(argv[0]) + 1);

  dolfin::init(petscargc, petscargv);

  while ((c = getopt_long(argc, argv, "hlpv:d:V", long_options, &option_index))!=-1)
  {
    switch (c)
    {
      case 'h':
        command_line_options["help"] = "";
        break;

      case 'l':
        command_line_options["log"] = "";
        break;

      case 'p':
        command_line_options["petsc-info"] = "";
        break;

      case 'v':
        command_line_options["verbose"] = optarg;
        break;

      case 'd':
        command_line_options["dolfin-verbose"] = optarg;
        break;

      case 'V':
        command_line_options["version"] = "";
        break;

      default:
        log(ERROR, "ERROR: Unrecognized option.");
        usage(argv[0]);
        std::exit(-1);
    }
  }

  if(command_line_options.count("verbose") == 0)                     // verbose
  {
    verbosity = WARNING;
  }
  else
  {
    verbosity = parse_verbosity(command_line_options["verbose"]);
  }

  set_log_level(verbosity);

  if(command_line_options.count("log"))                              // redirect std::err/std::out
  {                                                                  // FIXME: should be done by passing an ofstream to Logger
    std::ostringstream debug_file, err_file;                         // but this wouldn't redirect SPuD which automatically uses
    debug_file << "terraferma.log";                                  // std::cout
    err_file << "terraferma.err";
    int proc = dolfin::MPI::rank(MPI_COMM_WORLD);
    debug_file << "-" << proc;
    err_file << "-" << proc;

    if(std::freopen(debug_file.str().c_str(), "w", stdout) == NULL)
    {
      tf_err("Failed to redirect stdio for debugging.", "std::freopen failed.");
    }

    if(std::freopen(err_file.str().c_str(), "w", stderr) == NULL)
    {
      tf_err("Failed to redirect stderr for debugging.", "std::freopen failed.");
    }
  }

  if(command_line_options.count("help"))                             // help
  {
    usage(argv[0]);
    std::exit(-1);
  }

  if(command_line_options.count("version"))                          // version
  {
    print_version();
    std::exit(-1);
  }

  if(command_line_options.count("dolfin-verbose") == 0)              // dolfin-verbose
  {
    dolfinverbosity = verbosity;
  }
  else
  {
    dolfinverbosity = parse_dolfin_verbosity(command_line_options["verbose"]);
  }

  dolfin::set_log_active(true);
  dolfin::set_log_level(dolfinverbosity);

  if(command_line_options.count("petsc-info"))                       // petsc-info
  {
    PetscErrorCode perr;                                             // petsc error code
    perr = PetscInfoAllow(PETSC_TRUE, PETSC_NULL);
    petsc_err(perr);
  }

  if(argc > optind + 1)                                              // find the options file name
  {
    command_line_options["tfml"] = argv[optind + 1];
  }
  else if(argc == optind + 1)
  {
    command_line_options["tfml"] = argv[optind];
  }

  Spud::load_options(command_line_options["tfml"]);                  // and load it
  if(!Spud::have_option("/io/output_base_name"))                     // check its been loaded
  {
    tf_err("The input options file appears invalid.", "Failed to find required option /io/output_base_name in %s.", command_line_options["tfml"].c_str());
  }

  if(verbosity <= INFO)
  {
    std::stringstream s; s.str("");
    s << "Command line:" << std::endl;                       // print the command line
    for(int i=0;i<argc; i++)
    {
      s << argv[i] << " ";
    }
    s << std::endl;
    log(INFO, s.str());
    print_environment();                                             // print the environment variables

    Spud::print_options();                                           // print the options tree
  }

  return;
}

//*******************************************************************|************************************************************//
// return the git sha
//*******************************************************************|************************************************************//
std::string buckettools::githash()
{
  return std::string(__GIT_SHA__);
}

//*******************************************************************|************************************************************//
// return the TerraFERMA version
//*******************************************************************|************************************************************//
std::string buckettools::tfversion()
{
  return std::string(__TERRAFERMA_VERSION__);
}

//*******************************************************************|************************************************************//
// return the compile time
//*******************************************************************|************************************************************//
std::string buckettools::compiletime()
{
  std::stringstream s; s.str("");
  s << __DATE__<<" "<<__TIME__;
  return s.str();
}

