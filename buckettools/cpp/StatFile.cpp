
#include "StatFile.h"
#include "confdefs.h"
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

void StatFile::header_constants()
{
  std::string buffer;
  char * cbuffer;
  time_t rawtime;
  
  (*this).constant_tag("HGChangesetId", "string", __HG_ID__);
  buffer  = __DATE__;
  buffer += " ";
  buffer += __TIME__;
  (*this).constant_tag("CompileTime", "string", buffer);
  cbuffer = ctime(&rawtime);
//   if((*cbuffer+sizeof(&cbuffer)-1)=='\n')
//   {
//     (*cbuffer+sizeof(&cbuffer)-1) = '\0';
//   }
//   buffer = static_cast<std::string>(cbuffer);
  (*this).constant_tag("StartTime", "string", cbuffer);
  cbuffer = getenv("HOSTNAME");
  if(cbuffer==NULL)
  {
    
  }
//   buffer = static_cast<std::string>(cbuffer);
//   (*this).constant_tag("HostName", "string", cbuffer);
}

void StatFile::constant_tag(const std::string name, 
                            const std::string type, 
                            const std::string value)
{
  file_ << "<constant name=\""<<name
        <<"\" type=\""<<type
        <<"\" value=\""<<value<<"\" />" 
        << std::endl;
}


