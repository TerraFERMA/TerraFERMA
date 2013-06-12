
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
// print the recommended usage to std::err
//*******************************************************************|************************************************************//
void buckettools::usage(char *cmd)
{
  print_version(std::cerr);
  std::cerr << std::endl << std::endl << "Usage: " << cmd << " [options ...] [simulation-file]" << std::endl
      << std::endl << "Options:" << std::endl
      <<" -h, --help" << std::endl << "\tHelp! Prints this message."  << std::endl
      <<" -l, --log" << std::endl << "\tCreate log (redirects stdout) and error (redirects stderr) file for each process." << std::endl
      <<" -p, --petsc-info" << std::endl << "\tPETSc verbose output to stdout." << std::endl
      <<" -v <level>, --verbose" << std::endl << "\tVerbose output to stdout, default level dolfin::WARNING (30)." << std::endl
      <<" -V, --version" << std::endl << "\tVersion information ." << std::endl;
  return;
}

//*******************************************************************|************************************************************//
// print the version number to an std::ostream
//*******************************************************************|************************************************************//
void buckettools::print_version(std::ostream& stream)
{
  stream << "GitHash: " <<__GIT_SHA__ << std::endl
      <<"CompileTime: "<<__DATE__<<" "<<__TIME__<<std::endl
      ;
  stream.flush();
  return;
}

//*******************************************************************|************************************************************//
// print relevant environment variables to an std::ostream
//*******************************************************************|************************************************************//
void buckettools::print_environment(std::ostream& stream)
{

  const char *relevant_variables[]={"PETSC_OPTIONS"};

  int no_relevant_variables=1;
  int variables_found=0;

  stream<<"Relevant environment variables:" << std::endl;
  for(int i=0; i<no_relevant_variables; i++) {
    char *env=getenv(relevant_variables[i]);
    if(env!=NULL) {
      stream << relevant_variables[i]<<" = " << env <<std::endl;
      variables_found++;
    }
  };
  if (!variables_found)
  {
    stream<<" none" << std::endl;
  }
  stream.flush();
  return;
}


//*******************************************************************|************************************************************//
// parse the command line argument for verbosity to an integer
//*******************************************************************|************************************************************//
int buckettools::parse_verbosity(const std::string &verbosity)
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
    {"help",       no_argument,       0, 'h'},
    {"log",        no_argument,       0, 'l'},
    {"petsc-info", no_argument,       0, 'p'},
    {"verbose",    optional_argument, 0, 'v'},
    {"version",    no_argument,       0, 'V'},
    {0,            0,                 0, 0}                             // terminated with an array of zeros
  };
  int option_index = 0;
  int verbosity;
  int c;
  std::map< std::string, std::string > command_line_options;

  int petscargc;                                                    // dummy argc and argv to stop our command line arguments being
  char** petscargv;                                                 // interpretted by petsc usage.  FIXME: ugly hack!
  petscargc = 1;
  petscargv = new char*[1];
  petscargv[0] = new char[strlen(argv[0]) + 1];
  strncpy(petscargv[0], argv[0], strlen(argv[0]) + 1);

  dolfin::init(petscargc, petscargv);

  while ((c = getopt_long(argc, argv, "hlpv::V", long_options, &option_index))!=-1){
    switch (c){
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
      command_line_options["verbose"] = (optarg == NULL) ? "WARNING" : optarg;
      break;

    case 'V':
      command_line_options["version"] = "";
      break;

    default:
      dolfin::error("getopt returned unrecognized character code.");

    }
  }

  if(command_line_options.count("help"))                             // help
  {
    usage(argv[0]);
    std::exit(-1);
  }

  if(command_line_options.count("version"))                          // version
  {
    print_version(std::cerr);
    std::exit(-1);
  }

  if(command_line_options.count("verbose") == 0)                     // verbose
  {
    verbosity = dolfin::WARNING;
  }
  else
  {
    verbosity = parse_verbosity(command_line_options["verbose"]);
  }

  dolfin::set_log_active(true);
  dolfin::set_log_level(verbosity);

  if(command_line_options.count("petsc-info"))                       // petsc-info
  {
    PetscErrorCode perr;                                             // petsc error code
    perr = PetscInfoAllow(PETSC_TRUE, PETSC_NULL); CHKERRV(perr);
  }

  if(command_line_options.count("log"))                              // std::err/std::out
  {
    std::ostringstream debug_file, err_file;
    debug_file << "bucket.log";
    err_file << "bucket.err";
    int proc = dolfin::MPI::process_number();
    debug_file << "-" << proc;
    err_file << "-" << proc;

    if(std::freopen(debug_file.str().c_str(), "w", stdout) == NULL)
    {
      dolfin::error("Failed to redirect stdio for debugging");
    }

    if(std::freopen(err_file.str().c_str(), "w", stderr) == NULL)
    {
      dolfin::error("Failed to redirect stderr for debugging");
    }
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
    dolfin::error("Failed to find output_base_name after loading options file.");
  }

  if(verbosity <= dolfin::INFO)
  {
    std::cout << "Command line:" << std::endl;                       // print the command line
    for(int i=0;i<argc; i++)
    {
      std::cout<<argv[i]<<" ";
    }
    std::cout << std::endl;
    print_environment(std::cout);                                    // print the environment variables

    Spud::print_options();                                           // print the options tree to std::cout
  }

  return;
}


