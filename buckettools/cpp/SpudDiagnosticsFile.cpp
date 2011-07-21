
#include "SpudDiagnosticsFile.h"
#include "Bucket.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <spud>

using namespace buckettools;

SpudDiagnosticsFile::SpudDiagnosticsFile(const std::string name) : StatFile(name)
{
  // Do nothing... all handled by StatFile constructor
}

SpudDiagnosticsFile::~SpudDiagnosticsFile()
{
  // Do nothing... all handled by StatFile destructor
}

bool SpudDiagnosticsFile::include_in_file()
{

}


