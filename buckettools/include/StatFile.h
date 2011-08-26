
#ifndef __STAT_FILE_H
#define __STAT_FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include <dolfin.h>

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // StatFile class:
  //
  // This is a base class, which provides simple functionality, for writing various pieces of output to file in a parsable
  // format (using an xml header).  Actual files for specific purposes (diagnostics, detectors etc.) should be defined in
  // classes derived from this base class.
  //*****************************************************************|************************************************************//
  class StatFile
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    StatFile();                                                      // default constructor
    
    StatFile(const std::string name);                                // specific constructor
    
    virtual ~StatFile();                                             // default destructor
    
  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // file name

    std::ofstream file_;                                             // file stream

    //***************************************************************|***********************************************************//
    // Header writing functions
    //***************************************************************|***********************************************************//

    void header_constants_();                                        // write the header entries for constant values that do not 
                                                                     // appear later again in the file
    
    void header_timestep_(uint &column);                             // write the header for a dynamic simulation (including columns 
                                                                     // for the elapsed time, timestep etc.)
    
    void constant_tag_(const std::string &name,                      // write a header tag for a constant value
                       const std::string &type, 
                       const std::string &value);
    
    void tag_(const std::string &name,                               // write a header tag for a temporal value that does not belong
              const uint &column,                                    // to a system or have components
              const std::string &statistic)
    { tag_(name, column, statistic, "", 0); }

    void tag_(const std::string &name,                               // write a header tag for a temporal value that does not have
              const uint &column,                                    // components
              const std::string &statistic,
              const std::string &system)
    { tag_(name, column, statistic, system, 0); }
    
    void tag_(const std::string &name,                               // write a header tag for a temporal value (most generic)
              const uint &column,
              const std::string &statistic,
              const std::string &system,
              const uint &components);

    //***************************************************************|***********************************************************//
    // Data writing functions
    //***************************************************************|***********************************************************//

    void data_timestep_(const uint   &timestep,                      // write the data for timestepping for a dynamic simulation
                        const double &elapsedtime, 
                        const double &dt);
    
  };
  
}
#endif
