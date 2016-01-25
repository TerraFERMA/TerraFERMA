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


#ifndef __FUNCTIONALBUCKET_H
#define __FUNCTIONALBUCKET_H

#include "BoostTypes.h"
#include "PointDetectors.h"
#include "ReferencePoint.h"
#include <dolfin.h>

namespace buckettools
{


  //*****************************************************************|************************************************************//
  // FunctionalBucket class:
  //
  // The FunctionalBucket class describes a system functional.
  //*****************************************************************|************************************************************//
  class FunctionalBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    FunctionalBucket();                                                // default constructor

    FunctionalBucket(SystemBucket* system);                            // specific constructor
    
    virtual ~FunctionalBucket();                                       // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void attach_form_coeffs();

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    double value(const bool& force=false);                           // calculate and return the value of the functional

    double oldvalue() const
    { return oldvalue_; }

    double change();

    void update();

    void resetcalculated();

    const std::string name() const                                   // return a constant string giving the function name
    { return name_; }

    const std::string uflsymbol() const                              // return a constant string giving the ufl symbol 
    { return uflsymbol_; }                                           // for this function

    SystemBucket* system()                                           // return a pointer to the parent system
    { return system_; }

    const SystemBucket* system() const                               // return a constant pointer to the parent system
    { return system_; }

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    void output();                                                   // output mesh functions

    virtual const std::string str(const int indent=0) const;         // return an indented string describing the contents of this
                                                                     // functional

    virtual const bool include_in_statistics() const;                // return a boolean indicating if this function is included in 
                                                                     // diagnostic output

    virtual const bool include_in_steadystate() const;               // return a boolean indicating if this function is included in 
                                                                     // steady state output

    virtual const bool output_cellfunction() const;                  // return a boolean indicating if the functional is to 
                                                                     // be output as a cell function
    
    virtual const bool output_facetfunction() const;                 // return a boolean indicating if the functional is to 
                                                                     // be output as a facet function

  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // available to this class and its derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string name_;                                               // the functional name

    std::string uflsymbol_;                                          // the functional ufl symbol

    Form_ptr form_;                                                  // the form of the functional

    SystemBucket* system_;                                           // a pointer to the parent system

    bool calculated_;                                                // boolean indicating if the functional has been calculated

    double value_, oldvalue_;                                        // value and previous value

    dolfin::CellFunction<double> *cellfunction_;                     // cell function for output

    dolfin::FacetFunction<double> *facetfunction_;                   // facet function for output

  };

  typedef std::shared_ptr< FunctionalBucket > FunctionalBucket_ptr;    // define a (boost shared) pointer to the function bucket class type

}
#endif
