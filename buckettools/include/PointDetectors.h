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


#ifndef __POINT_DETECTORS_H
#define __POINT_DETECTORS_H

#include <dolfin.h>
#include "GenericDetectors.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // PointDetectors class:
  //
  // PointDetectors is a derived class of GenericDetectors that implements a detector at a single point
  //*****************************************************************|************************************************************//
  class PointDetectors : public GenericDetectors
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    PointDetectors();                                                // generic constructor
    
    PointDetectors(const Array_double_ptr point, 
                                    const std::string &name);        // specific constructor (dolfin array)
    
    PointDetectors(const std::vector<double> &point, 
                                        const std::string &name);    // specific constructor (std vector)
    
    ~PointDetectors();                                               // default destructor
    
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:
    
    //***************************************************************|***********************************************************//
    // Initialization
    //***************************************************************|***********************************************************//

    void init_(const Array_double_ptr point);                        // initialize the point detector at a point defined by an array
    
    void init_(const std::vector<double> &point);                    // initialize the point detector at a point defined by a vector
    
  };
  
  typedef boost::shared_ptr< PointDetectors > PointDetectors_ptr;    // define a (boost shared) pointer for this class type
  
}

#endif
