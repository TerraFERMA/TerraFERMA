
#ifndef __SIGINTEVENTHANDLER_H
#define __SIGINTEVENTHANDLER_H

#include <dolfin.h>
#include <signal.h>
#include "EventHandler.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SigIntEventHandler class:
  //
  // A derived class that instructs the signal handler how to react to sigints
  //*****************************************************************|************************************************************//

  class SigIntEventHandler : public EventHandler
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    SigIntEventHandler();                                            // default constructor

    ~SigIntEventHandler();                                           // default destructor

    virtual int handle_signal(int signum);                           // hook method

    sig_atomic_t received()                                          // accessor
    { return received_; }

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    sig_atomic_t received_;                                          // whether a sigint has been received or not

  };

  typedef boost::shared_ptr< SigIntEventHandler > 
                                             SigIntEventHandler_ptr; // define a boost shared ptr type for the class
}
#endif

