
#include "StatFile.h"
#include "builddefs.h"
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace buckettools;

StatFile::StatFile()
{
  // Do nothing
}

StatFile::StatFile(const std::string name) : name_(name)
{
  // Open the file member
  file_.open((char*)name.c_str());
}

StatFile::~StatFile()
{
  // maybe close the _file member here?
  if (file_.is_open())
  {
    file_.close();
  }
}

void StatFile::header_constants_()
{
  std::string buffer;
  char * cbuffer;
  time_t rawtime;
  
  constant_tag_("HGChangesetId", "string", __HG_ID__);
  
  buffer  = __DATE__;
  buffer += " ";
  buffer += __TIME__;
  constant_tag_("CompileTime", "string", buffer);
  
//   cbuffer = ctime(&rawtime);
//   if((*cbuffer+sizeof(&cbuffer)-1)=='\n')
//   {
//     (*cbuffer+sizeof(&cbuffer)-1) = '\0';
//   }
//   buffer = static_cast<std::string>(cbuffer);
//   constant_tag_("StartTime", "string", cbuffer);
  
//   cbuffer = getenv("HOSTNAME");
//   if(cbuffer==NULL)
//   {
//     
//   }
//   buffer = static_cast<std::string>(cbuffer);
//   constant_tag_("HostName", "string", cbuffer);
}

void StatFile::header_timestep_(uint &column)
{
  
  tag_("timestep", column++, "value");
  tag_("ElapsedTime", column++, "value");
  tag_("dt", column++, "value");
  
}

void StatFile::constant_tag_(const std::string &name, 
                             const std::string &type, 
                             const std::string &value)
{
  
  file_ << "<constant name=\"" << name
        << "\" type=\"" << type
        << "\" value=\"" << value << "\" />" 
        << std::endl << std::flush;

}

void StatFile::tag_(const std::string &name,
                    const uint &column,
                    const std::string &statistic,
                    const std::string &system,
                    const uint &components)
{
  
  file_ << "<field column=\"" << column
        << "\" name=\"" << name
        << "\" statistic=\"" << statistic << "\"";
  if(!system.empty())
  {
    file_ << " system=\"" << system << "\"";
  }
  if (components > 0)
  {
    file_ << " components=\"" << components << "\"";
  }
  file_ << " />" << std::endl << std::flush;

  
}


