
#include <dolfin.h>
#include "petscsnes.h"
#include "BucketPETScBase.h"

using namespace buckettools;

PetscErrorCode buckettools::FormFunction(SNES snes, Vec x, Vec f, void* ctx)
{
  SNESCtx      *snesctx    = (SNESCtx *) ctx;
  Function_ptr iteratedfunction = (*snesctx).iteratedfunction;
  std::vector<BoundaryCondition_ptr>& bcs = (*snesctx).bcs;
  Vec_ptr px(&x, dolfin::NoDeleter());
  Vec_ptr pf(&f, dolfin::NoDeleter());
  dolfin::PETScVector rhs(pf), iteratedvec(px);
  bool reset_tensor = false;

  std::cout << "In FormFunction" << std::endl;

  (*iteratedfunction).vector() = iteratedvec;

  dolfin::assemble(rhs, (*(*snesctx).linear), reset_tensor);
  for(uint i = 0; i < bcs.size(); ++i)
  {
    (*bcs[i]).apply(rhs, iteratedvec);
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode buckettools::FormJacobian(SNES snes, Vec x, Mat *A, Mat *B, MatStructure* flag, void* ctx)
{
  SNESCtx *snesctx = (SNESCtx *) ctx;
  Function_ptr iteratedfunction = (*snesctx).iteratedfunction;
  std::vector<BoundaryCondition_ptr>& bcs = (*snesctx).bcs;
  Vec_ptr px(&x, dolfin::NoDeleter());
  Mat_ptr pA(A, dolfin::NoDeleter());
  Mat_ptr pB(B, dolfin::NoDeleter());
  dolfin::PETScVector iteratedvec(px);
  dolfin::PETScMatrix matrix(pA), matrixpc(pB);
  bool reset_tensor = false;

  std::cout << "In FormJacobian" << std::endl;

  (*iteratedfunction).vector() = iteratedvec;

  dolfin::assemble(matrix, (*(*snesctx).bilinear), reset_tensor);
  for(uint i = 0; i < bcs.size(); ++i)
  {
    (*bcs[i]).apply(matrix);
  }
  if ((*snesctx).bilinearpc)
  {
    dolfin::assemble(matrixpc, (*(*snesctx).bilinearpc), reset_tensor);
    for(uint i = 0; i < bcs.size(); ++i)
    {
      (*bcs[i]).apply(matrixpc);
    }
  }

  *flag = SAME_NONZERO_PATTERN;

  PetscFunctionReturn(0);
}


