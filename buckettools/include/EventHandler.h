
#ifndef __EVENTHANDLER_H
#define __EVENTHANDLER_H

#include <dolfin.h>
#include <signal.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // EventHandler class:
  //
  // A base class that instructs the signal handler how to react to different signals
  //*****************************************************************|************************************************************//

  class EventHandler
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    virtual int handle_signal(int signum) = 0;                       // pure virtual hook for signals

    virtual sig_atomic_t received() = 0;                             // pure virtual received check for signals

  };

  typedef boost::shared_ptr< EventHandler > EventHandler_ptr;        // define a boost shared ptr type for the class

}
#endif

