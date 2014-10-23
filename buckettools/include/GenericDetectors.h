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


#ifndef __GENERICDETECTORS_H
#define __GENERICDETECTORS_H

#include <dolfin.h>
#include "BoostTypes.h"

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // GenericDetectors class:
  //
  // GenericDetectors is a generic class that contains base information for point and python detectors
  //*****************************************************************|************************************************************//
  class GenericDetectors
  {
  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:                                                            // accessible to everyone

    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    GenericDetectors();                                              // default constructor
    
    GenericDetectors(const uint &number_detectors,                   // specific constructor - no. detectors, coordinate dimension,
                     const uint &meshdim,                            // and name
                     const std::string &name);
    
    virtual ~GenericDetectors();                                     // destructor
    
    //***************************************************************|***********************************************************//
    // Detector evaluation
    //***************************************************************|***********************************************************//

    void eval(std::vector< Array_double_ptr > &values,               // base eval implementation, take values of a function at
              const dolfin::GenericFunction &function,               // detector positions and returns values
              Mesh_ptr mesh);                                   

    void eval_ownership(Mesh_ptr mesh);                              // evaluate and store the cell and detector ownership of 
                                                                     // detectors on a mesh

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string name() const                                   // return the name
    { return name_; }

    const uint dim() const                                           // return the coordinate dimensions
    { return meshdim_; } 

    const uint size() const                                          // return the number of detectors in this set
    { return number_detectors_; }

    const std::vector<int> cell_ids(Mesh_ptr mesh);                  // return the cell ids

    const std::vector<int> detector_ids(Mesh_ptr mesh);              // return the detector ids (owned by this process)

    //***************************************************************|***********************************************************//
    // Data array access
    //***************************************************************|***********************************************************//

    std::vector< Array_double_ptr >::iterator begin();               // return an iterator to the beginning of the positions

    std::vector< Array_double_ptr >::const_iterator begin() const;   // return a constant iterator to the beginning of the positions

    std::vector< Array_double_ptr >::iterator end();                 // return an iterator to the end of the positions

    std::vector< Array_double_ptr >::const_iterator end() const;     // return a constant iterator to the end of the positions
    
    //***************************************************************|***********************************************************//
    // Output
    //***************************************************************|***********************************************************//

    virtual const std::string str() const;                           // return a string giving the detector locations
    
  //*****************************************************************|***********************************************************//
  // Protected functions
  //*****************************************************************|***********************************************************//

  protected:                                                         // only available to this class and derived classes
    
    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string                     name_;                           // detectors set name
    std::size_t                     number_detectors_,               // number of detectors
                                    meshdim_;                        // coordinate dimensions
    std::vector< Array_double_ptr > positions_;                      // vector of arrays giving locations of detectors

    std::map< Mesh_ptr, std::vector< int > > cell_ids_;              // the cell ids for a particular mesh - not initialized until eval is called
    std::map< Mesh_ptr, std::vector< int > > detector_ids_;          // the detectors ids that this process owns - not initialized until eval is called
    
    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void clean_();                                                   // empty the data maps
    
  };
  
  typedef std::shared_ptr< GenericDetectors > GenericDetectors_ptr;  // define a (boost shared) pointer to the class type
  
}

#endif
