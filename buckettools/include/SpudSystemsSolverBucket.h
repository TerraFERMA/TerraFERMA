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


#ifndef __SPUD_SYSTEMS_SOLVER_BUCKET_H
#define __SPUD_SYSTEMS_SOLVER_BUCKET_H

#include "BoostTypes.h"
#include "SystemsSolverBucket.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudSystemsSolverBucket class:
  //
  //*****************************************************************|************************************************************//
  class SpudSystemsSolverBucket : public SystemsSolverBucket
  {
  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // available to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudSystemsSolverBucket();                                       // default constructor

    SpudSystemsSolverBucket(const std::string &optionpath, 
                            const int &solve_location, 
                            Bucket* bucket, 
                            SystemsSolverBucket* systemssolver=NULL);// specific constructor (with parent bucket and solver)
    
    virtual ~SpudSystemsSolverBucket();                                  // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill in the data in the base solver bucket

    void initialize_diagnostics();                                   // initialize any diagnostic output from this systems solver

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a constant string containing the optionpath for the solver bucket
    { return optionpath_; }

    const std::string str(int indent=0) const;                       // return an indented string describing the contents of this systems solver

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible by this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the optionpath for the solver bucket

    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill the base data of the systems solver bucket

    void fill_solvers_();                                            // fill the child solvers

    
    //***************************************************************|***********************************************************//
    // Output functions (continued)
    //***************************************************************|***********************************************************//

    void checkpoint_options_();                                      // checkpoint the options system for the spudsystembucket

  };
 
  typedef std::shared_ptr< SpudSystemsSolverBucket > SpudSystemsSolverBucket_ptr;// define a (boost shared) pointer to a spud systems solver bucket class

}
#endif

