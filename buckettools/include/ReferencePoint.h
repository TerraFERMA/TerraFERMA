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


#ifndef __REFERENCEPOINT_H
#define __REFERENCEPOINT_H

#include <dolfin.h>
#include "BoostTypes.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // ReferencePoint class:
  //
  // ReferencePoint implements a method for setting reference points on functions.
  //*****************************************************************|************************************************************//
  class ReferencePoint : public dolfin::DirichletBC
  {
  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    ReferencePoint(const std::vector<double> &coord, 
                   const FunctionSpace_ptr functionspace,
                   const GenericFunction_ptr value);
    
    ~ReferencePoint();                                              // default destructor

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    static SubDomain_ptr subdomain_(const std::vector<double> &coord,
                                    const FunctionSpace_ptr functionspace);

  };

  typedef std::shared_ptr< ReferencePoint > ReferencePoint_ptr;    // define a (boost shared) pointer for this class type

  //*****************************************************************|************************************************************//
  // RerencePointSubDomain class:
  //
  // An overloaded dolfin subdomain class to describe the side of an internally generated dolfin mesh 
  //*****************************************************************|************************************************************//
  class ReferencePointSubDomain : public dolfin::SubDomain
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accesible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    ReferencePointSubDomain(const std::vector<double> &point, 
                            const double &tolerance=1.e-10);          // optional constructor

    ~ReferencePointSubDomain();                                                         // default destructor

    //***************************************************************|***********************************************************//
    // Overloaded base class functions
    //***************************************************************|***********************************************************//
    
    bool inside(const dolfin::Array<double>& x, bool on_boundary)    // return a boolean, true if point x is inside the subdomain, false otherwise
                                                              const; 

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible in this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::vector<double> point_;                                      // point we want to be near

  };
}

#endif
