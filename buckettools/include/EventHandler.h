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

