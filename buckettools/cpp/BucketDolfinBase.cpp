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
#include "Logger.h"
#include <dolfin.h>
#include <fstream>

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

std::string buckettools::xml_filename(const std::string& basename)
{
  // This routine tests whether various extensions exist and returns a filename
  // for an xml file with an appropriate extension

  std::string filename;
  std::ifstream file;
  std::stringstream filenamestream;

  std::vector<std::string> extensions = {"", ".gz", ".xml", ".xml.gz"};
  std::vector<std::string>::const_iterator ext_it;
  for (ext_it = extensions.begin(); ext_it != extensions.end(); ++ext_it)
  {
    filenamestream.str(""); filenamestream << basename << *ext_it;
    file.open(filenamestream.str().c_str(), std::ifstream::in);
    if (file)
    {
      file.close();
      filename = filenamestream.str();
      break;
    }
  }

  if (ext_it == extensions.end())
  {
    tf_err("Could not find requested xml file.",
           "%s, %s.gz, %s.xml or %s.xml.gz not found.", 
           basename.c_str(), basename.c_str(), basename.c_str(), basename.c_str());
  }

  return filename;
}

