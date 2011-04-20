
#ifndef __STAT_FILE_H
#define __STAT_FILE_H

#include <cstdio>
#include <fstream>
#include <string>

namespace buckettools
{

  class StatFile
  {
  public:
    
    /// Default Constructor
    StatFile();
    
    /// Constructor - takes the _file _name
    StatFile(const std::string name);
    
    /// Destructor - closes the file
    ~StatFile();
    
  protected:
    
    /// name of the file
    std::string _name;
    /// the file itself
    std::ofstream _file;
  };
  
}
#endif