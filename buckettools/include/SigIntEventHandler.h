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

    void handle_signal(int signum);                                   // hook method

    sig_atomic_t received()                                          // accessor
    { return received_; }

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:

    sig_atomic_t received_;                                          // whether a sigint has been received or not

  };

  typedef std::shared_ptr< SigIntEventHandler > 
                                             SigIntEventHandler_ptr; // define a boost shared ptr type for the class
}
#endif

