
#include <dolfin.h>
#include "petscsnes.h"
#include "BucketPETScBase.h"
#include "SolverBucket.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// define the petsc snes callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::FormFunction(SNES snes, Vec x, Vec f, 
                                                          void* ctx)
{

  dolfin::log(dolfin::INFO, "In FormFunction");

  SNESCtx *snesctx = (SNESCtx *)ctx;                                 // cast the snes context

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver

  if ((*solver).monitor_norms())
  {
    PetscErrorCode perr;                                             // petsc error code
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(1): 2-norm x = %f", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(1): inf-norm x = %f", norm);

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(1): 2-norm f = %f", norm);

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(1): inf-norm f = %f", norm);
  }

  Function_ptr iteratedfunction = (*snesctx).iteratedfunction;       // collect the iterated system bucket function
  std::vector<BoundaryCondition_ptr>& bcs = (*snesctx).bcs;          // get the vector of bcs
  std::vector<ReferencePoints_ptr>& points = (*snesctx).points;      // get the vector of reference points
  Vec_ptr px(&x, dolfin::NoDeleter());                               // convert the iterated snes vector
  Vec_ptr pf(&f, dolfin::NoDeleter());                               // convert the rhs snes vector
  dolfin::PETScVector rhs(pf), iteratedvec(px);
  bool reset_tensor = false;                                         // never reset the tensor

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*(*snesctx).bucket).update_nonlinear();                           // update nonlinear coefficients

  dolfin::assemble(rhs, (*(*snesctx).linear), reset_tensor);         // assemble the rhs from the context linear form
  for(uint i = 0; i < bcs.size(); ++i)                               // loop over the bcs
  {
    (*bcs[i]).apply(rhs, iteratedvec);                               // FIXME: will break symmetry?
  }
  
  for(uint i = 0; i < points.size(); ++i)                            // loop over the reference points
  {
    (*points[i]).apply(rhs, iteratedvec);                            // FIXME: will break symmetry?
  }
  
  if ((*solver).monitor_norms())
  {
    PetscErrorCode perr;                                             // petsc error code
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(2): 2-norm x = %f", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(2): inf-norm x = %f", norm);

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(2): 2-norm f = %f", norm);

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormFunction(2): inf-norm f = %f", norm);
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc snes callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::FormJacobian(SNES snes, Vec x, Mat *A, 
                                         Mat *B, MatStructure* flag, 
                                         void* ctx)
{
  dolfin::log(dolfin::INFO, "In FormJacobian");

  SNESCtx *snesctx = (SNESCtx *)ctx;                                 // cast the snes context

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver


  if ((*solver).monitor_norms())
  {
    PetscErrorCode perr;                                             // petsc error code
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): 2-norm x = %f", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): inf-norm x = %f", norm);

    perr = MatNorm(*A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): Frobenius norm A = %f", norm);

    perr = MatNorm(*A,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): inf-norm A = %f", norm);

    perr = MatNorm(*B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): Frobenius norm B = %f", norm);

    perr = MatNorm(*B,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(1): inf-norm B = %f", norm);
  }

  Function_ptr iteratedfunction = (*snesctx).iteratedfunction;       // collect the iterated system bucket function
  std::vector<BoundaryCondition_ptr>& bcs = (*snesctx).bcs;          // get the vector of bcs
  std::vector<ReferencePoints_ptr>& points = (*snesctx).points;      // get the vector of reference points
  Vec_ptr px(&x, dolfin::NoDeleter());                               // convert the iterated snes vector
  dolfin::PETScVector iteratedvec(px);
  Mat_ptr pA(A, dolfin::NoDeleter());                                // convert the snes matrix
  Mat_ptr pB(B, dolfin::NoDeleter());                                // convert the snes matrix pc
  dolfin::PETScMatrix matrix(pA), matrixpc(pB);
  bool reset_tensor = false;                                         // never reset the tensor

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*(*snesctx).bucket).update_nonlinear();                           // update nonlinear coefficients

  dolfin::assemble(matrix, (*(*snesctx).bilinear), reset_tensor);    // assemble the matrix from the context bilinear form
  for(uint i = 0; i < bcs.size(); ++i)                               // loop over the bcs
  {
    (*bcs[i]).apply(matrix);                                         // FIXME: will break symmetry?
  }
  for(uint i = 0; i < points.size(); ++i)                            // loop over the reference points
  {
    (*points[i]).apply(matrix);                                      // FIXME: will break symmetry?
  }
  if ((*snesctx).ident_zeros)
  {
    matrix.ident_zeros();
  }

  if ((*snesctx).bilinearpc)                                         // do we have a different bilinear pc form associated?
  {
    dolfin::assemble(matrixpc, (*(*snesctx).bilinearpc),             // assemble the matrix pc from the context bilinear pc form
                                                      reset_tensor);
    for(uint i = 0; i < bcs.size(); ++i)                             // loop over the bcs
    {
      (*bcs[i]).apply(matrixpc);                                     // FIXME: will break symmetry
    }
    for(uint i = 0; i < points.size(); ++i)                          // loop over the points
    {
      (*points[i]).apply(matrixpc);                                  // FIXME: will break symmetry
    }
    if ((*snesctx).ident_zeros_pc)
    {
      matrixpc.ident_zeros();
    }
  }

  *flag = SAME_NONZERO_PATTERN;                                      // both matrices are assumed to have the same sparsity

  if ((*solver).monitor_norms())
  {
    PetscErrorCode perr;                                             // petsc error code
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): 2-norm x = %f", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): inf-norm x = %f", norm);

    perr = MatNorm(*A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): Frobenius norm A = %f", norm);

    perr = MatNorm(*A,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): inf-norm A = %f", norm);

    perr = MatNorm(*B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): Frobenius norm B = %f", norm);

    perr = MatNorm(*B,NORM_INFINITY,&norm); CHKERRQ(perr);
    dolfin::log(dolfin::get_log_level(), "FormJacobian(2): inf-norm B = %f", norm);
  }

  PetscFunctionReturn(0);
}


