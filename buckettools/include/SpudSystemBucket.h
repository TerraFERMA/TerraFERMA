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


#ifndef __SPUD_SYSTEM_H
#define __SPUD_SYSTEM_H

#include "BoostTypes.h"
#include "SystemBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudSystemBucket class:
  //
  // The SpudSystemBucket class is a derived class of the system that populates the
  // data structures within a system using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudSystemBucket : public SystemBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudSystemBucket(const std::string &optionpath, Bucket* bucket); // default constructor

    ~SpudSystemBucket();                                             // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill the system assuming a buckettools spud schema

    void allocate_coeff_function();                                  // allocate the coefficient functions

    void initialize_fields_and_coefficient_expressions();            // initialize the expressions that fields and coefficients use

    void initialize_coefficient_functions();                         // initialize the coefficients functions

    void initialize_solvers();                                       // initialize the matrices related to forms as well as the
                                                                     // petsc objects

    void initialize_diagnostics() const;                             // initialize any diagnostic output from this system

    std::vector< GenericFunction_ptr > collect_vis_functions() const;// output the diagnostics on this system

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a string containing the optionpath for the system
    { return optionpath_; }
    
    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str(int indent=0) const;                       // return an indented string describing the contents of the
                                                                     // system

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible in this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the system optionpath

    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill in the base information about the system

    void fill_systemfunction_();                                     // fill in the system function information

    void fill_fields_();                                             // fill in the data about the system fields (subfunctions)

    void fill_bcs_();                                                // fill in the data about the system bcs

    void fill_coeffs_();                                             // fill in the coefficient information

    void fill_solvers_();                                            // fill in the solver bucket information

    void fill_functionals_();                                        // fill in the functional bucket information

  };
 
  typedef std::shared_ptr< SpudSystemBucket > SpudSystemBucket_ptr;// define a (boost shared) pointer type for the spud system

}
#endif
