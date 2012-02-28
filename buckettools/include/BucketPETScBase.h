
#ifndef __BUCKETPETSC_BASE_H
#define __BUCKETPETSC_BASE_H

#include "BoostTypes.h"
#include "petscsnes.h"
#include "Bucket.h"

namespace buckettools
{

  //*****************************************************************|************************************************************//
  // A collection of structures and callback functions used in the bucket by petsc
  //*****************************************************************|************************************************************//

  typedef struct {                                                   // a structure used to pass bucket data into snes callback functions
    Form_ptr linear;                                                 // linear form
    Form_ptr bilinear;                                               // bilinear form
    bool ident_zeros;                                                // the matrix has zeros rows which should be identified
    Form_ptr bilinearpc;                                             // bilinear pc (may be null if not used)
    bool ident_zeros_pc;                                             // the pc has zero rows which should be identified
    std::vector<BoundaryCondition_ptr> bcs;                          // bcs
    Function_ptr iteratedfunction;                                   // work function
    Bucket *bucket;                                                  // pointer to bucket
    SolverBucket *solver;                                            // pointer to solver
  } SNESCtx;

  PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx);   // petsc snes callback function to form the residual

  PetscErrorCode FormJacobian(SNES snes, Vec x, Mat *A, Mat *B,      // petsc snes callback function to form the jacobian
                                   MatStructure* flag, void* ctx);

}

#endif
