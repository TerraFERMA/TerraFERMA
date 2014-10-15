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


#ifndef __BUCKETPETSC_BASE_H
#define __BUCKETPETSC_BASE_H

#include "BoostTypes.h"
#include "petscsnes.h"
#include "Bucket.h"
#include "ConvergenceFile.h"
#include "KSPConvergenceFile.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // A collection of structures and callback functions used in the bucket by petsc
  //*****************************************************************|************************************************************//

  typedef struct {                                                   // a structure used to pass bucket data into snes callback functions
    SolverBucket *solver;                                            // pointer to solver
  } SNESCtx;

  PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx);   // petsc snes callback function to form the residual

  PetscErrorCode FormJacobian(SNES snes, Vec x, Mat *A, Mat *B,      // petsc snes callback function to form the jacobian
                                   MatStructure* flag, void* ctx);

  typedef struct {                                                   // a structure used to pass bucket data into monitor functions
    SolverBucket *solver;                                            // pointer to solver
  } CustomMonitorCtx;

  PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt its,          // petsc snes callback function to output a 
                                      PetscReal norm, void* mctx);   // convergence file

  PetscErrorCode KSPCustomMonitor(KSP ksp, int it,                   // petsc ksp callback function to output a 
                                      PetscReal rnorm, void* mctx);  // convergence file

  PetscErrorCode KSPNullSpaceMonitor(KSP ksp, int it,                // petsc ksp callback function to test a
                                     PetscReal rnorm, void* mctx);   // null space

  PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu);

  void petsc_failure(PetscErrorCode perr,
                     const std::string &filename,
                     const int &line,
                     const std::string &dirname,
                     const std::string &petsc_function);

  #define petsc_fail(perr) do {petsc_failure(perr, __FILE__, __LINE__, __SDIR__, PETSC_FUNCTION_NAME);} while(0)

  void petsc_error(PetscErrorCode perr,
                   const std::string &filename,
                   const int &line,
                   const std::string &dirname,
                   const std::string &petsc_function);

  #define petsc_err(perr)  do {petsc_error(perr, __FILE__, __LINE__, __SDIR__, PETSC_FUNCTION_NAME);} while(0)

}

#endif
