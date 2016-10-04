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


#ifndef __SPUD_FUNCTIONALBUCKET_H
#define __SPUD_FUNCTIONALBUCKET_H

#include "BoostTypes.h"
#include "FunctionalBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudFunctionalBucket class:
  //
  // The SpudFunctionalBucket class is a derived class of the functional that populates the
  // data structures within a functional using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudFunctionalBucket : public FunctionalBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudFunctionalBucket(const std::string &optionpath, 
                                               SystemBucket* system);// specific constructor
    
    ~SpudFunctionalBucket();                                           // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill this functional

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a string containing the optionpath of this function
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str(const int indent=0) const;                 // return an indented string describing the contents of the functional bucket

    const bool include_in_statistics() const;                        // return a boolean indicating if this function is to 
                                                                     // be included in diagnostic output
    
    const bool include_in_steadystate() const;                       // return a boolean indicating if this function is to 
                                                                     // be included in steadystate output
    
    const bool output_cellfunction() const;                          // return a boolean indicating if the functional is to 
                                                                     // be output as a cell function
    
    const bool output_facetfunction() const;                         // return a boolean indicating if the functional is to 
                                                                     // be output as a facet function
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only available this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the optionpath of this function

  };
 
  typedef std::shared_ptr< SpudFunctionalBucket >                    // define a (boost shared) pointer for this function class type
                                            SpudFunctionalBucket_ptr;

}
#endif
