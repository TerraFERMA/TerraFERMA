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


#include "SignalHandler.h"
#include "EventHandler.h"
#include <dolfin.h>
#include <string>
#include <signal.h>

using namespace buckettools;

SignalHandler *SignalHandler::instance_ = NULL;                      // initialize the global static class variables
EventHandler_ptr SignalHandler::signal_handlers_[NSIG];

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SignalHandler::SignalHandler()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// register a handler
//*******************************************************************|************************************************************//
EventHandler_ptr SignalHandler::register_handler(int signum, EventHandler_ptr eh)
{

  EventHandler_ptr old_eh = SignalHandler::signal_handlers_[signum]; // copy the old_eh from the signum slot

  SignalHandler::signal_handlers_[signum] = eh;                      // store eh in the signum slot in the table
 
  struct sigaction sa;                                               // register the dispatcher to handle this signum
  sa.sa_handler = SignalHandler::dispatcher;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sigaction (signum, &sa, 0);

  return old_eh;
}

//*******************************************************************|************************************************************//
// return a handler
//*******************************************************************|************************************************************//
EventHandler_ptr SignalHandler::return_handler(int signum)
{

  return SignalHandler::signal_handlers_[signum];                    // return a pointer to the corresponding entry in the table
                                                                     // (may be NULL)
}


//*******************************************************************|************************************************************//
// dispatch the command requested for this signal
//*******************************************************************|************************************************************//
void SignalHandler::dispatcher(int signum)
{
  if (SignalHandler::signal_handlers_[signum] != 0)                  // sanity check
  {
    (*SignalHandler::signal_handlers_[signum]).handle_signal (signum);
  }
}

//*******************************************************************|************************************************************//
// return the signal handler
//*******************************************************************|************************************************************//
SignalHandler* SignalHandler::instance()
{
  if(!instance_)
  {
    instance_ = new SignalHandler;
  }
  return SignalHandler::instance_;
}

