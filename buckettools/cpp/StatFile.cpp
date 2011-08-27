
#include "StatFile.h"
#include "builddefs.h"
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
StatFile::StatFile()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
StatFile::StatFile(const std::string &name) : name_(name)
{
  file_.open((char*)name.c_str());                                   // open the file_ member
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
StatFile::~StatFile()
{
  if (file_.is_open())
  {
    file_.close();                                                   // close the file_ member
  }
}

//*******************************************************************|************************************************************//
// write lines of the xml header for constants that do not vary throughout a simulation
//*******************************************************************|************************************************************//
void StatFile::header_constants_()
{
  std::string buffer;
  char * cbuffer;
  time_t rawtime;
  
  constant_tag_("HGChangesetId", "string", __HG_ID__);               // the mercurial changeset id
  
  buffer  = __DATE__;
  buffer += " ";
  buffer += __TIME__;
  constant_tag_("CompileTime", "string", buffer);                    // the comilation time
  
//   cbuffer = ctime(&rawtime);
//   if((*cbuffer+sizeof(&cbuffer)-1)=='\n')
//   {
//     (*cbuffer+sizeof(&cbuffer)-1) = '\0';
//   }
//   buffer = static_cast<std::string>(cbuffer);
//   constant_tag_("StartTime", "string", cbuffer);                      // the simulation start time
  
//   cbuffer = getenv("HOSTNAME");
//   if(cbuffer==NULL)
//   {
//     
//   }
//   buffer = static_cast<std::string>(cbuffer);
//   constant_tag_("HostName", "string", cbuffer);                       // the hostname
}

//*******************************************************************|************************************************************//
// write lines of the xml header for values relating to timestepping
//*******************************************************************|************************************************************//
void StatFile::header_timestep_(uint &column)
{
  
  tag_("timestep", column++, "value");                               // the timestep count (i.e. number of timesteps taken)
  tag_("ElapsedTime", column++, "value");                            // the current time
  tag_("dt", column++, "value");                                     // the actual timestep
  
}

//*******************************************************************|************************************************************//
// write an xml tag for a constant that does not vary throughout a simulation
//*******************************************************************|************************************************************//
void StatFile::constant_tag_(const std::string &name, 
                             const std::string &type, 
                             const std::string &value)
{
  
  file_ << "<constant name=\"" << name
        << "\" type=\"" << type
        << "\" value=\"" << value << "\" />" 
        << std::endl << std::flush;

}

//*******************************************************************|************************************************************//
// write an xml tag for a variable in a simulation 
//*******************************************************************|************************************************************//
void StatFile::tag_(const std::string &name,
                    const uint &column,
                    const std::string &statistic,
                    const std::string &system,
                    const uint &components)
{
  
  file_ << "<field column=\"" << column
        << "\" name=\"" << name
        << "\" statistic=\"" << statistic << "\"";
  if(!system.empty())                                                // is this part of a system?
  {
    file_ << " system=\"" << system << "\"";
  }
  if (components > 0)                                                // does it have subcomponents (i.e. is it rank>0)? 
  {
    file_ << " components=\"" << components << "\"";
  }
  file_ << " />" << std::endl << std::flush;

  
}

//*******************************************************************|************************************************************//
// write data to the file for values relating to timestepping
//*******************************************************************|************************************************************//
void StatFile::data_timestep_(const uint   &timestep,
                              const double &elapsedtime, 
                              const double &dt)
{
  
  file_.setf(std::ios::scientific);
  file_.precision(10);
  
  file_ << timestep << " ";  
  file_ << elapsedtime << " ";
  file_ << dt << " ";
  
  file_.unsetf(std::ios::scientific);
  
}

