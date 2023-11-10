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


#ifndef __GENERICSOLVERBUCKET_H
#define __GENERICSOLVERBUCKET_H

#include "BoostTypes.h"

namespace buckettools
{
  
  enum solve_location { SOLVE_START, SOLVE_TIMELOOP, SOLVE_DIAGNOSTICS, SOLVE_NEVER };

  //*****************************************************************|************************************************************//
  // GenericSolverBucket class:
  //
  //*****************************************************************|************************************************************//
  class GenericSolverBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Functions used to run the model
    //***************************************************************|***********************************************************//

    virtual bool solve() = 0;                                        // run the nonlinear solver described by this class

    virtual double residual_norm() = 0;                              // return the norm of the residual (which will be reassembled)

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return a string containing the solver name
    { return name_; }

    const bool solved() const                                        // return a boolean indicating if this system solver has been solved
    { return *solved_; }                                             // for or not

    virtual void set_current_systemssolver(const std::string systemssolvername) {}

    virtual void reset_current_systemssolver() {}

    virtual void resetcalculated() = 0;                              // update this solver at the end of a timestep

    //***************************************************************|***********************************************************//
    // Output
    //***************************************************************|***********************************************************//

    virtual const bool include_in_visualization() const = 0;

    virtual const bool include_in_statistics() const = 0;

    virtual const bool include_in_steadystate() const = 0;

    virtual const bool include_in_detectors() const = 0;

    virtual const std::string str(int indent=0) const = 0;           // return an indented string describing the contents of this solver

  protected:                                                         // availble to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // solver name

    bool_ptr solved_;                                                // indicate if the system has been solved this timestep


  };

  typedef std::shared_ptr< GenericSolverBucket > GenericSolverBucket_ptr;// define a boost shared ptr type for the class

}
#endif

