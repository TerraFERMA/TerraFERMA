
#ifndef __DIAGNOSTICS_FILE_H
#define __DIAGNOSTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "StatFile.h"

namespace buckettools
{

  class DiagnosticsFile : public StatFile
  {
  public:
    
    DiagnosticsFile(const std::string name);
    
    void write_header();

  };
  
}
#endif