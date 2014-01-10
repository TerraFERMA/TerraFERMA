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


#ifndef __BUCKETDOLFIN_BASE_H
#define __BUCKETDOLFIN_BASE_H

#include <dolfin.h>

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // A collection of classes and subroutines used in the bucket by dolfin
  //*****************************************************************|************************************************************//

  //*****************************************************************|************************************************************//
  // Side class:
  //
  // An overloaded dolfin subdomain class to describe the side of an internally generated dolfin mesh 
  //*****************************************************************|************************************************************//
  class Side : public dolfin::SubDomain
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accesible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    Side(const uint &component, const double &side);                 // optional constructor

    ~Side();                                                         // default destructor

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

    uint component_;                                                 // component of n-dimensional space

    double side_;                                                    // location of the side in that n-dimensional space

  };

}

bool abslessthan(const double &elem1, const double &elem2);

#endif

