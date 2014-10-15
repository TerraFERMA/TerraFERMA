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


#include "Logger.h"
#include "Usage.h"
#include "SignalHandler.h"

using namespace buckettools;

static boost::scoped_array<char> va_buffer_(0);
static unsigned int va_buffer_size_= 0;

// Buffer allocation
void allocate_va_buffer(std::string msg)
{
  // vsnprintf requires a char pointer of fixed size so we
  // need to allocate the va_buffer here. We allocate twice the size of
  // the format string and at least 1000.
  unsigned int new_size = std::max(static_cast<unsigned int>(2*msg.size()),
                                   static_cast<unsigned int>(1000));
  //static_cast<unsigned int>(DOLFIN_LINELENGTH));
  if (new_size > va_buffer_size_)
  {
    va_buffer_.reset(new char[new_size]);
    va_buffer_size_ = new_size;
  }
}

// Macro for parsing arguments
#define va_read(va_buffer, msg) \
  allocate_va_buffer(msg); \
  va_list va_ptr; \
  va_start(va_ptr, msg); \
  vsnprintf(va_buffer, va_buffer_size_, msg.c_str(), va_ptr); \
  va_end(va_ptr);

Logger* Logger::instance_ = NULL;                                    // initialize the global static class variables

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
Logger::Logger() : logstream_(&std::cout), errstream_(&std::cerr),
                   loglevel_(WARNING)
{
                                                                     // do nothing
}

Logger::Logger(const Logger &logger)
{
  loglevel_ = logger.loglevel_;
  logstream_ = logger.logstream_;
  errstream_ = logger.errstream_;
}

//*******************************************************************|************************************************************//
// set log stream
//*******************************************************************|************************************************************//
void Logger::set_log_stream(std::ostream& logstream)
{
  logstream_ = &logstream;
}

//*******************************************************************|************************************************************//
// set err stream
//*******************************************************************|************************************************************//
void Logger::set_err_stream(std::ostream& errstream)
{
  errstream_ = &errstream;
}

//*******************************************************************|************************************************************//
// set log level
//*******************************************************************|************************************************************//
void Logger::set_log_level(int &loglevel)
{
  loglevel_ = loglevel;
}

//*******************************************************************|************************************************************//
// throw a sigint
//*******************************************************************|************************************************************//
void Logger::failure(const std::string &filename, 
                     const int &line,
                     const std::string &errstr,
                     const std::string &reason) const
{
  std::stringstream s;
  s << "*** ERROR: terminating at the end of the timestep."
    << std::endl 
    << description(filename, line, errstr, reason)
    << std::endl;

  write(ERROR, s.str());

  (*SignalHandler::instance()).dispatcher(SIGINT);
}

//*******************************************************************|************************************************************//
// throw a runtime error
//*******************************************************************|************************************************************//
void Logger::error(const std::string &filename, 
                   const int &line,
                   const std::string &errstr,
                   const std::string &reason) const
{
  std::stringstream s;
  s << "*** ERROR: terminating immediately!"
    << std::endl
    << description(filename, line, errstr, reason)
    << std::endl;

  write(ERROR, s.str());

  throw std::runtime_error("std::runtime_error thrown.");
}

//*******************************************************************|************************************************************//
// print a warning
//*******************************************************************|************************************************************//
void Logger::warning(const std::string &filename, 
                     const int &line,
                     const std::string &errstr,
                     const std::string &reason) const
{
  std::stringstream s;
  s << std::endl
    << "*** WARNING:" 
    << std::endl
    << description(filename, line, errstr, reason)
    << std::endl;

  write(WARNING, s.str());
}

//*******************************************************************|************************************************************//
// return a description of the problem/error/warning
//*******************************************************************|************************************************************//
std::string Logger::description(const std::string &filename, 
                                const int &line,
                                const std::string &errstr,
                                const std::string &reason) const
{
  std::stringstream s;
  s << "-------------------------------------------------------------------------"
    << std::endl
    << "*** " << "Problem:  " << errstr << std::endl
    << "*** " << "Reason:   " << reason << std::endl
    << "*** " << "Where:    " << filename << ":" << line << std::endl
    << "*** " << std::endl
    << "*** " << "GitHash:  " << githash() << std::endl
    << "*** " << "Compiled: " << compiletime() << std::endl
    << "*** " << std::endl
    << "-------------------------------------------------------------------------"
    << std::endl;

  return s.str();
}

//*******************************************************************|************************************************************//
// write to the log or the error output (depending on the loglevel)
//*******************************************************************|************************************************************//
void Logger::write(int loglevel, std::string msg) const
{
  if (loglevel<loglevel_)
  {
    return;
  }

  if (loglevel>=WARNING)
  {
    *errstream_ << msg << std::endl;
  }
  else
  {
    *logstream_ << msg << std::endl;
  }
}

//*******************************************************************|************************************************************//
// return the signal handler
//*******************************************************************|************************************************************//
Logger* Logger::instance()
{
  if(!instance_)
  {
    instance_ = new Logger;
  }
  return instance_;
}

void buckettools::set_log_level(int &loglevel)
{
  (*Logger::instance()).set_log_level(loglevel);
}

void buckettools::log(int loglevel, std::string msg, ...)
{
  va_read(va_buffer_.get(), msg);
  (*Logger::instance()).write(loglevel, va_buffer_.get());
}

void buckettools::failure(const std::string &filename, 
                          const int &line,
                          const std::string &errstr,
                          const std::string &reason, ...)
{
  va_read(va_buffer_.get(), reason);
  (*Logger::instance()).failure(filename, line, errstr, va_buffer_.get());
}
                            
void buckettools::error(const std::string &filename, 
                        const int &line,
                        const std::string &errstr,
                        const std::string &reason, ...)
{
  va_read(va_buffer_.get(), reason);
  (*Logger::instance()).error(filename, line, errstr, va_buffer_.get());
}
                            
void buckettools::warning(const std::string &filename, 
                          const int &line,
                          const std::string &errstr,
                          const std::string &reason, ...)
{
  va_read(va_buffer_.get(), reason);
  (*Logger::instance()).warning(filename, line, errstr, va_buffer_.get());
}
                            
