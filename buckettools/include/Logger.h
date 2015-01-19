// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __LOGGER_H
#define __LOGGER_H

#include <ostream>

namespace buckettools
{

  enum LogLevel
  {
    ERROR     = 40, // things that go boom
    WARNING   = 30, // things that may go boom later
    INFO      = 20, // information of general interest
    DBG       = 10  // sundry
  };
  
  //*****************************************************************|************************************************************//
  // Logger class:
  //
  // A class that logs.
  //*****************************************************************|************************************************************//
  class Logger
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    static Logger* instance();                                       // entry point - no public constructor in a singleton class

    void set_log_stream(std::ostream& logstream);

    void set_err_stream(std::ostream& errstream);

    void set_log_level(int &loglevel);

    void failure(const std::string &filename, 
                 const int &line,
                 const std::string &errstr,
                 const std::string &reason) const;

    void error(const std::string &filename, 
               const int &line,
               const std::string &errstr,
               const std::string &reason) const;

    void warning(const std::string &filename, 
                 const int &line,
                 const std::string &errstr,
                 const std::string &reason) const;

    std::string description(const std::string &filename, 
                            const int &line, 
                            const std::string &errstr, 
                            const std::string &reason) const;

    void write(int loglevel, std::string msg) const;


  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to members of this class

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    Logger();                                                        // a private constructor so a singleton

    Logger(const Logger& logger);

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    static Logger *instance_;                                        // singleton pointer

    std::ostream *logstream_, *errstream_;                            // log and error streams

    int loglevel_;

  };

  void log(int loglevel, std::string msg, ...);

  void set_log_level(int &loglevel);

  void failure(const std::string &filename, 
               const int &line,
               const std::string &errstr,
               const std::string &reason, ...);

  #define tf_fail(errstr, reason, ...) do {failure(__FILE__, __LINE__, errstr, reason,##__VA_ARGS__);} while(0)

  void error(const std::string &filename, 
             const int &line,
             const std::string &errstr,
             const std::string &reason, ...);

  #define tf_err(errstr, reason, ...) do {error(__FILE__, __LINE__, errstr, reason,##__VA_ARGS__);} while(0)

  void warning(const std::string &filename, 
               const int &line,
               const std::string &errstr,
               const std::string &reason, ...);

  #define tf_warn(errstr, reason, ...) do {warning(__FILE__, __LINE__, errstr, reason,##__VA_ARGS__);} while(0)

}

#endif

