
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

}

#endif
