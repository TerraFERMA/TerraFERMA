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


#include "PointDetectors.h"
#include "GenericDetectors.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
PointDetectors::PointDetectors() : GenericDetectors()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PointDetectors::PointDetectors(const Array_double_ptr point, 
                                          const std::string &name) : 
                           GenericDetectors(1, (*point).size(), name)
{
  init_(point);                                                      // initialize the detector 
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
PointDetectors::PointDetectors(const std::vector<double> &point, 
                                          const std::string &name) : 
                          GenericDetectors(1, point.size(), name)
{
  init_(point);                                                      // initialize the detector
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
PointDetectors::~PointDetectors()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// intialize a point detector from a (boost shared) pointer to a dolfin array
//*******************************************************************|************************************************************//
void PointDetectors::init_(const Array_double_ptr point)
{
  if(!positions_.empty())
  {
    dolfin::error("In PointDetectors::init_ intializing already initialized detectors.");
  }
  
  positions_.push_back(point);
}

//*******************************************************************|************************************************************//
// intialize a point detector from a std vector
//*******************************************************************|************************************************************//
void PointDetectors::init_(const std::vector<double> &point)
{
  Array_double_ptr arraypoint(new dolfin::Array<double>(point.size()));
  for (uint i = 0; i<point.size(); i++)
  {
    (*arraypoint)[i] = point[i];
  }
  
  init_(arraypoint);
}

