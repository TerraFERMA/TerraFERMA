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


#include "BoostTypes.h"
#include "FunctionalBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include "BucketDolfinBase.h"
#include "DolfinPETScBase.h"
#include "BucketPETScBase.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
FunctionalBucket::FunctionalBucket() : value_(0.0), oldvalue_(0.0), calculated_(false)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
FunctionalBucket::FunctionalBucket(SystemBucket* system) : system_(system), value_(0.0), oldvalue_(0.0), calculated_(false)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
FunctionalBucket::~FunctionalBucket()
{
  if (cellfunction_)
  {
    delete cellfunction_;
  }
  if (facetfunction_)
  {
    delete facetfunction_;
  }
}

//*******************************************************************|************************************************************//
// attach the coefficients this functional requests using the parent bucket data maps
//*******************************************************************|************************************************************//
void FunctionalBucket::attach_form_coeffs()
{
  (*(*system_).bucket()).attach_coeffs(form_);
}

//*******************************************************************|************************************************************//
// return the value of the functional (calculating it if necessary)
//*******************************************************************|************************************************************//
double FunctionalBucket::value(const bool& force)
{
  if(!calculated_ || force)
  {

    if (output_cellfunction())
    {
      if (!cellfunction_)
      {
        cellfunction_ = new dolfin::CellFunction<double>((*system()).mesh());
      }
      (*cellfunction_).set_all(0.0);
    }

    if (output_facetfunction())
    {
      if (!facetfunction_)
      {
        facetfunction_ = new dolfin::FacetFunction<double>((*system()).mesh());
      }
      (*facetfunction_).set_all(0.0);
    }

    dolfin::Assembler assembler;
    dolfin::Scalar value;
    assembler.assemble(value, *form_, cellfunction_, facetfunction_);
    value_ = value.get_scalar_value();
    if (!force)                                                      // if this has been forced it is outside the normal scope
    {                                                                // of us calculating the functional value for output
      calculated_ = true;                                            // so don't update the calculated flag or it'll throw off
    }                                                                // our output assumptions
  }
  return value_;
}

//*******************************************************************|************************************************************//
// return the values of the functional per cell
//*******************************************************************|************************************************************//
dolfin::CellFunction<double> FunctionalBucket::cellfunction(const bool& force) const
{
  // set up new functions so this can be const
  dolfin::CellFunction<double> l_cellfunction((*system()).mesh());

  // only calculate the cell function if we haven't calculated it yet or we're not
  // outputing the cell function (in which case it won't be included in our normal
  // calculation of the functional anyway) or if we're forcing it
  if (!force && calculated_ && output_cellfunction())
  {
    l_cellfunction = *cellfunction_;
  }
  else
  {
    l_cellfunction.set_all(0.0);
    dolfin::FacetFunction<double>* l_facetfunction = NULL;
  
    dolfin::Assembler assembler;
    dolfin::Scalar value;
    assembler.assemble(value, *form_, &l_cellfunction, l_facetfunction);

  }

  return l_cellfunction;
}

//*******************************************************************|************************************************************//
// return the values of the functional per facet
//*******************************************************************|************************************************************//
dolfin::FacetFunction<double> FunctionalBucket::facetfunction(const bool& force) const
{
  // set up new functions so this can be const
  dolfin::FacetFunction<double> l_facetfunction((*system()).mesh());

  // only calculate the cell function if we haven't calculated it yet or we're not
  // outputing the cell function (in which case it won't be included in our normal
  // calculation of the functional anyway) or if we're forcing it
  if (!force && calculated_ && output_cellfunction())
  {
    l_facetfunction = *facetfunction_;
  }
  else
  {
    l_facetfunction.set_all(0.0);
    dolfin::CellFunction<double>* l_cellfunction = NULL;
  
    dolfin::Assembler assembler;
    dolfin::Scalar value;
    assembler.assemble(value, *form_, l_cellfunction, &l_facetfunction);

  }

  return l_facetfunction;
}

//*******************************************************************|************************************************************//
// return the change in the value of the functional over a timestep
//*******************************************************************|************************************************************//
double FunctionalBucket::change()
{
  double fvalue = value();
  return std::abs(fvalue-oldvalue())/std::max(std::abs(fvalue), DOLFIN_EPS);
}

//*******************************************************************|************************************************************//
// update the timelevels of the functional
//*******************************************************************|************************************************************//
void FunctionalBucket::update()
{
  if (include_in_steadystate())
  {
    double fvalue = value();                                         // check that the value has been calculated this timestep
  }
  oldvalue_ = value_;
}

//*******************************************************************|************************************************************//
// output the functional (calculating it if necessary)
//*******************************************************************|************************************************************//
void FunctionalBucket::output()
{
  if (output_cellfunction() || output_facetfunction())
  {
    if (!calculated_)
    {
      double fvalue = value();
    }

    if (output_cellfunction())
    {
      assert(cellfunction_);
      
      std::stringstream buffer;
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_"
                             << name() << "_"
                             << (*(*system()).bucket()).visualization_count() << "_cellfunction.xml";
      dolfin::File cellfunction_file(buffer.str());
      cellfunction_file << *cellfunction_;
    }

    if (output_facetfunction())
    {
      assert(facetfunction_);
      
      std::stringstream buffer;
      buffer.str(""); buffer << (*(*system()).bucket()).output_basename() << "_" 
                             << (*system()).name() << "_"
                             << name() << "_" 
                             << (*(*system()).bucket()).visualization_count() << "_facetfunction.xml";
      dolfin::File facetfunction_file(buffer.str());
      facetfunction_file << *facetfunction_;
    }
  }

}

//*******************************************************************|************************************************************//
// reset calculated flag
//*******************************************************************|************************************************************//
void FunctionalBucket::resetcalculated()
{
  calculated_ = false;
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the functional bucket
//*******************************************************************|************************************************************//
const std::string FunctionalBucket::str(const int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionalBucket " << name() << std::endl;
  return s.str();
}

//*******************************************************************|************************************************************//
// include this function in diagnostic output
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionalBucket::include_in_statistics() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_statistics.");
  return false;
}

//*******************************************************************|************************************************************//
// include this function in steadystate output and checking
// this is a virtual function and should be implemented in the derived options class
//*******************************************************************|************************************************************//
const bool FunctionalBucket::include_in_steadystate() const
{
  tf_err("Failed to find virtual function.", "Need a virtual include_in_steadystate.");
  return false;
}

//*******************************************************************|************************************************************//
// output the cell integrals of the functional as a cell function
//*******************************************************************|************************************************************//
const bool FunctionalBucket::output_cellfunction() const
{
  tf_err("Failed to find virtual function.", "Need a virtual output_cellfunction.");
  return false;
}

//*******************************************************************|************************************************************//
// output the cell integrals of the functional as a facet function
//*******************************************************************|************************************************************//
const bool FunctionalBucket::output_facetfunction() const
{
  tf_err("Failed to find virtual function.", "Need a virtual output_facetfunction.");
  return false;
}



