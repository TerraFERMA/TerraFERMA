
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

    static void dispatcher(int signum);                              // entry point adapter installled into sigaction

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

