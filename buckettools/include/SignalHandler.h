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


#ifndef __SIGNALHANDLER_H
#define __SIGNALHANDLER_H

#include <dolfin.h>
#include <signal.h>
#include "EventHandler.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SignalHandler class:
  //
  // A class that allows signals to be intercepted and custom functions written to control their behaviour.
  //*****************************************************************|************************************************************//
  class SignalHandler
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    static SignalHandler *instance();                                // entry point - no public constructor in a singleton class

    EventHandler_ptr register_handler (int signum,                   // register an event handler for signum and return a pointer to
                        EventHandler_ptr eh);                        // any existing event handler that was previously registered
                                                                     // for signum

    EventHandler_ptr return_handler(int signum);                     // return a pointer to the handler associated with signum

    static void dispatcher(int signum);                              // entry point adapter installed into sigaction

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible to members of this class

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SignalHandler();                                                 // a private constructor so a singleton

    static SignalHandler *instance_;                                 // singleton pointer

    static EventHandler_ptr signal_handlers_[NSIG];                  // table ot pointers to eventhandlers registered by applications
                                                                     // NSIG is the number of signals defined in signal.h
  };

}

#endif

