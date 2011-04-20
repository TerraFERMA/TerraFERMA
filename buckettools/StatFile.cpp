
#include "StatFile.h"
#include "hgid.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>

using namespace buckettools;

StatFile::StatFile()
{
  // Do nothing
}

StatFile::StatFile(const std::string name) : _name(name)
{
  // Open the file member
  _file.open((char*)name.c_str());
}

StatFile::~StatFile()
{
  // maybe close the _file member here?
  if (_file.is_open())
  {
    _file.close();
  }
}


