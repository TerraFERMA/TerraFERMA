
#ifndef __STAT_FILE_H
#define __STAT_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include <dolfin.h>

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
    virtual ~StatFile();
    
  protected:
    
    /// name of the file
    std::string name_;
    /// the file itself
    std::ofstream file_;


    void header_constants_();
    
    void header_timestep_(uint &column);
    
    void constant_tag_(const std::string &name, 
                       const std::string &type, 
                       const std::string &value);
    
    void tag_(const std::string &name,
              const uint &column,
              const std::string &statistic)
    { tag_(name, column, statistic, "", 0); }

    void tag_(const std::string &name,
              const uint &column,
              const std::string &statistic,
              const std::string &system)
    { tag_(name, column, statistic, system, 0); }
    
    void tag_(const std::string &name,
              const uint &column,
              const std::string &statistic,
              const std::string &system,
              const uint &components);
    
  };
  
}
#endif
