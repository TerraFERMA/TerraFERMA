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


#include "BucketDolfinBase.h"
#include <dolfin.h>

using namespace buckettools;

Side::Side(const uint &component, const double &side) : component_(component), side_(side)
{
  // Do nothing
}

Side::~Side()
{
  // Do nothing
}

bool Side::inside(const dolfin::Array<double>& x, bool on_boundary) const
{
  return (std::fabs(x[component_] - side_) < DOLFIN_EPS && on_boundary);
}

