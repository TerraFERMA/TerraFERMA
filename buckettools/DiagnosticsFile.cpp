
#include "DiagnosticsFile.h"
#include "hgid.h"
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
  _file << "<header>" << std::endl;
  _file << "<constant name=\"HGChangesetId\" type=\"string\" value=\""<<__HG_ID__<<"\" />" << std::endl;
  _file << "</header>" << std::endl;
}


