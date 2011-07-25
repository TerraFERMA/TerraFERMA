
#ifndef __BUCKETPETSC_BASE_H
#define __BUCKETPETSC_BASE_H

#include "BoostTypes.h"
#include "petscsnes.h"

namespace buckettools
{

  typedef struct {
    //pointers to forms and boundary conditions
    Form_ptr linear;
    Form_ptr bilinear;
    Form_ptr bilinearpc;
    std::vector<BoundaryCondition_ptr>& bcs;

    Function_ptr function;    // work function
    Function_ptr oldfunction; // previous time step

  } SNESCtx;

  PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx);

  PetscErrorCode FormJacobian(SNES snes, Vec x, Mat *A, Mat *B, MatStructure* flag, void* ctx);

}

#endif
