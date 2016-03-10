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


#include "PythonExpression.h"
#include "BoostTypes.h"
#include "SpudFunctionalBucket.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemFunctionalsWrapper.h"
#include "SystemExpressionsWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "RegionsExpression.h"
#include "SemiLagrangianExpression.h"
#include "PointDetectors.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudFunctionalBucket::SpudFunctionalBucket(const std::string &optionpath, 
                                            SystemBucket* system) : 
                                            optionpath_(optionpath), 
                                            FunctionalBucket(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudFunctionalBucket::~SpudFunctionalBucket()
{
}

//*******************************************************************|************************************************************//
// fill the function bucket data structures assuming the buckettools schema and that this is a field
//*******************************************************************|************************************************************//
void SpudFunctionalBucket::fill()
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // field or coefficient name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/ufl_symbol";           // ufl symbol
  serr = Spud::get_option(buffer.str(), uflsymbol_); 
  spud_err(buffer.str(), serr);

  form_ = ufc_fetch_functional((*system_).name(), name_,             // get a pointer to the functional form from the ufc
                               (*system_).mesh());
  (*form_).set_cell_domains((*system_).celldomains());
  (*form_).set_interior_facet_domains((*system_).facetdomains());
  (*form_).set_exterior_facet_domains((*system_).facetdomains());

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // functional because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
  uint ncoeff = (*form_).num_coefficients();                         // how many coefficients does this functional require?
  for (uint i = 0; i < ncoeff; i++)
  {
    std::string uflsymbol = (*form_).coefficient_name(i);            // what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
    if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))    // a base ufl symbol was only inserted into the parent bucket's
    {                                                                // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
      std::string baseuflsymbol =                                    // what is the base ufl symbol?
            (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol);   // have we already registered a functionspace for this base ufl
                                                                     // symbol?
      if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
      {                                                              // no...
        FunctionSpace_ptr coefficientspace;
        coefficientspace = 
                      ufc_fetch_coefficientspace_from_functional(    // take a pointer to the functionspace from the ufc
                                      (*system_).name(), name(), 
                                      baseuflsymbol, 
                                      (*system_).mesh());
        (*(*system_).bucket()).register_coefficientspace(            // and register it in the parent bucket's map
                                      coefficientspace, 
                                      baseuflsymbol);
      }
    }

  }

  cellfunction_ = NULL;
  facetfunction_ = NULL;

}

//*******************************************************************|************************************************************//
// return a string describing the contents of the functional bucket
//*******************************************************************|************************************************************//
const std::string SpudFunctionalBucket::str(const int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "FunctionalBucket " << name() << " (" << 
                                    optionpath() << ")" << std::endl;
  return s.str();
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in diagnostic output
//*******************************************************************|************************************************************//
const bool SpudFunctionalBucket::include_in_statistics() const
{
  return Spud::have_option(optionpath()+"/include_in_statistics");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if this function bucket should be included in steady state check and output
//*******************************************************************|************************************************************//
const bool SpudFunctionalBucket::include_in_steadystate() const
{
  return Spud::have_option(optionpath()+"/include_in_steady_state");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the named functional of this function bucket should have a cell function output
//*******************************************************************|************************************************************//
const bool SpudFunctionalBucket::output_cellfunction() const
{
  return Spud::have_option(optionpath()+"/output_cell_function");
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the named functional of this function bucket should have a facet function output
//*******************************************************************|************************************************************//
const bool SpudFunctionalBucket::output_facetfunction() const
{
  return Spud::have_option(optionpath()+"/output_facet_function");
}

