
#ifndef __SPUD_DIAGNOSTICS_FILE_H
#define __SPUD_DIAGNOSTICS_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include "StatFile.h"
#include "Bucket.h"

namespace buckettools
{

  class SpudDiagnosticsFile : public StatFile
  {
  public:
    
    SpudDiagnosticsFile(const std::string name);
 
    ~SpudDiagnosticsFile();
    
    void write_header(const Bucket &bucket, 
                      const bool &timestepping);

  private:

    void header_bucket_(const Bucket &bucket,
                        uint &column);

  };
  
}
#endif
