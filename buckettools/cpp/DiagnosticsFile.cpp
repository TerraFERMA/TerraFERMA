
#include "DiagnosticsFile.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>

using namespace buckettools;

DiagnosticsFile::DiagnosticsFile(const std::string name) : StatFile(name)
{
  // Do nothing... all handled by StatFile constructor
}

void DiagnosticsFile::write_header()
{
  file_ << "<header>" << std::endl;
  header_constants_();
  file_ << "</header>" << std::endl;
}


