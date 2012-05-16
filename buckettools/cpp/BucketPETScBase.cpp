
#include <dolfin.h>
#include "petscsnes.h"
#include "BucketPETScBase.h"
#include "BucketDolfinBase.h"
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
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

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

  Function_ptr iteratedfunction = (*system).iteratedfunction();      // collect the iterated system bucket function
  const std::vector< const dolfin::DirichletBC* >& bcs = 
                                         (*system).dirichletbcs();   // get the vector of bcs
  const std::vector<ReferencePoints_ptr>& points = (*system).points();// get the vector of reference points
  Vec_ptr px(&x, dolfin::NoDeleter());                               // convert the iterated snes vector
  Vec_ptr pf(&f, dolfin::NoDeleter());                               // convert the rhs snes vector
  dolfin::PETScVector rhs(pf), iteratedvec(px);
  bool reset_tensor = false;                                         // never reset the tensor

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::assemble(rhs, *(*solver).linear_form(), reset_tensor);     // assemble the rhs from the context linear form
  for(uint i = 0; i < bcs.size(); ++i)                               // loop over the bcs
  {
    (*bcs[i]).apply(rhs, iteratedvec);
  }
  
  for(uint i = 0; i < points.size(); ++i)                            // loop over the reference points
  {
    (*points[i]).apply(rhs, iteratedvec);
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

  PetscErrorCode perr;                                               // petsc error code

  SNESCtx *snesctx = (SNESCtx *)ctx;                                 // cast the snes context

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  if ((*solver).monitor_norms())
  {
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

  Function_ptr iteratedfunction = (*system).iteratedfunction();      // collect the iterated system bucket function
  const std::vector< const dolfin::DirichletBC* >& bcs = 
                                         (*system).dirichletbcs();   // get the vector of bcs
  const std::vector<ReferencePoints_ptr>& points = (*system).points();// get the vector of reference points
  Vec_ptr px(&x, dolfin::NoDeleter());                               // convert the iterated snes vector
  dolfin::PETScVector iteratedvec(px);
  Mat_ptr pA(A, dolfin::NoDeleter());                                // convert the snes matrix
  Mat_ptr pB(B, dolfin::NoDeleter());                                // convert the snes matrix pc
  dolfin::PETScMatrix matrix(pA), matrixpc(pB);

  PETScMatrix_ptr matrixbc = (*solver).matrixbc();

  bool reset_tensor = false;                                         // never reset the tensor

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::symmetric_assemble(matrix, *matrixbc, 
                             *(*solver).bilinear_form(), bcs,
                             NULL, NULL, NULL,
                             reset_tensor, false, false);            // assemble the matrix from the context bilinear form
  Assembler::add_zeros_diagonal(matrix);
  matrix.apply("add");
  Assembler::add_zeros_diagonal(*matrixbc);
  (*matrixbc).apply("add");
  for(uint i = 0; i < points.size(); ++i)                            // loop over the reference points
  {
    (*points[i]).apply(matrix);
  }
  if ((*solver).ident_zeros())
  {
    matrix.ident_zeros();
  }

  if ((*solver).bilinearpc_form())                                   // do we have a different bilinear pc form associated?
  {
    dolfin::symmetric_assemble(matrixpc, *matrixbc, 
                               (*(*solver).bilinearpc_form()), bcs,
                               NULL, NULL, NULL,                     // assemble the matrix pc from the context bilinear pc form
                               reset_tensor, false, false);
    Assembler::add_zeros_diagonal(matrixpc);
    matrixpc.apply("add");
    Assembler::add_zeros_diagonal(*matrixbc);
    (*matrixbc).apply("add");
    for(uint i = 0; i < points.size(); ++i)                          // loop over the points
    {
      (*points[i]).apply(matrixpc);
    }
    if ((*solver).ident_zeros_pc())
    {
      matrixpc.ident_zeros();
    }
  }

  *flag = SAME_NONZERO_PATTERN;                                      // both matrices are assumed to have the same sparsity

  if ((*solver).monitor_norms())
  {
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

//*******************************************************************|************************************************************//
// define the petsc snes monitor callback function that outputs a visualization file
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::SNESCustomMonitor(SNES snes, PetscInt its,
                                      PetscReal norm, void* mctx)
{
  dolfin::log(dolfin::INFO, "In SNESCustomMonitor");

  std::stringstream buffer;                                          // string buffer
  PetscErrorCode perr;                                               // petsc error code

  CustomMonitorCtx *snesctx = (CustomMonitorCtx *)mctx;

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  (*solver).iteration_count(its);                                    // set the iteration count

  Vec x;
  perr = SNESGetSolution(snes, &x);  CHKERRQ(perr);                  // get the solution vector from snes
  Vec_ptr px(&x, dolfin::NoDeleter());                               // and convert it into a dolfin PETScVector
  dolfin::PETScVector sol(px);

  *(*(*system).iteratedfunction()).vector() = sol;

  Vec rx;
  perr = SNESGetFunction(snes, &rx, PETSC_NULL, PETSC_NULL);         // get the residual function vector from snes
  CHKERRQ(perr);
  Vec_ptr prx(&rx, dolfin::NoDeleter());                             // and convert it into a dolfin PETScVector
  dolfin::PETScVector solresid(prx);

  *(*(*system).residualfunction()).vector() = solresid;

  Vec dx;
  perr = SNESGetSolutionUpdate(snes, &dx);  CHKERRQ(perr);           // get the solution update vector from snes
  Vec_ptr pdx(&dx, dolfin::NoDeleter());                             // and convert it into a dolfin PETScVector
  dolfin::PETScVector solupdate(pdx);

  *(*(*system).snesupdatefunction()).vector() = solupdate;


  if ((*solver).convergence_file())
  {

    (*(*solver).convergence_file()).write_data();

  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ksp monitor callback function that outputs a convergence file
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::KSPCustomMonitor(KSP ksp, int it,
                                      PetscReal rnorm, void* mctx)
{
  dolfin::log(dolfin::INFO, "In KSPCustomMonitor");

  std::stringstream buffer;                                          // string buffer
  PetscErrorCode perr;                                               // petsc error code

  CustomMonitorCtx *kspctx = (CustomMonitorCtx *)mctx;

  SolverBucket* solver = (*kspctx).solver;                           // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  Vec x;
  perr = KSPBuildSolution(ksp, PETSC_NULL, &x);  CHKERRQ(perr);      // get the solution vector from the ksp
  Vec_ptr px(&x, dolfin::NoDeleter());                               // and convert it into a dolfin PETScVector
  dolfin::PETScVector sol(px);

  *(*(*system).iteratedfunction()).vector() = sol;

  Vec rx;
  perr = KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &rx);  CHKERRQ(perr);// get the residual vector from the ksp
  Vec_ptr prx(&rx, dolfin::NoDeleter());                             // and convert it into a dolfin PETScVector
  dolfin::PETScVector solresid(prx);

  *(*(*system).residualfunction()).vector() = solresid;

  if ((*solver).type()=="SNES")
  {
    PetscInt iter;
    perr = SNESGetIterationNumber((*solver).snes(), &iter); CHKERRQ(perr);
    (*solver).iteration_count(iter);
  }

  if ((*solver).ksp_convergence_file())
  {

    (*(*solver).ksp_convergence_file()).write_data(it);

  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ksp monitor callback function that tests the null space
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::KSPNullSpaceMonitor(KSP ksp, int it,
                                                PetscReal rnorm, void* mctx)
{
  if (it==0)
  {
    dolfin::log(dolfin::INFO, "In KSPNullSpaceMonitor");

    PetscErrorCode perr;                                             // petsc error code

    MatNullSpace SP;

    perr = KSPGetNullSpace(ksp, &SP); CHKERRQ(perr);

    if (SP)
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      PetscBool isNull;
      #else
      PetscTruth isNull;
      #endif

      Mat Amat, Pmat;
      MatStructure flag;

      perr = KSPGetOperators(ksp, &Amat, &Pmat, &flag); CHKERRQ(perr);

      perr = MatNullSpaceTest(SP, Amat, &isNull); CHKERRQ(perr);

      if (!isNull)
      {
        dolfin::error("Provided null space is not a null space of the matrix.");
      }
    }

  }

  PetscFunctionReturn(0);
}



