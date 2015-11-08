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


#include <dolfin.h>
#include "petscsnes.h"
#include "BucketPETScBase.h"
#include "Bucket.h"
#include "SystemBucket.h"
#include "SolverBucket.h"
#include "Logger.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// define the petsc snes callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::FormFunction(SNES snes, Vec x, Vec f, 
                                                          void* ctx)
{
  log(INFO, "In FormFunction");

  SNESCtx *snesctx = (SNESCtx *)ctx;                                 // cast the snes context

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  PetscErrorCode perr;                                               // petsc error code
  if ((*solver).monitor_norms())
  {
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(1): 2-norm x = %g", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(1): inf-norm x = %g", norm);

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(1): 2-norm f = %g", norm);

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(1): inf-norm f = %g", norm);
  }

  Function_ptr iteratedfunction = (*system).iteratedfunction();      // collect the iterated system bucket function
  const std::vector< const dolfin::DirichletBC* >& bcs = 
                                         (*system).bcs();   // get the vector of bcs
  dolfin::PETScVector rhs(f), iteratedvec(x);

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::Assembler assembler;
  assembler.assemble(rhs, *(*solver).linear_form());
  for(uint i = 0; i < bcs.size(); ++i)                               // loop over the bcs
  {
    (*bcs[i]).apply(rhs, (*(*iteratedfunction).vector()));
  }
  
  MatNullSpace sp = (*solver).nullspace();
  if (sp)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNullSpaceRemove(sp, rhs.vec(), PETSC_NULL);
    #else
    perr = MatNullSpaceRemove(sp, rhs.vec());
    #endif
    CHKERRQ(perr);
  }
  
  if ((*solver).monitor_norms())
  {
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(2): 2-norm x = %g", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(2): inf-norm x = %g", norm);

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(2): 2-norm f = %g", norm);

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormFunction(2): inf-norm f = %g", norm);
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc snes callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
PetscErrorCode buckettools::FormJacobian(SNES snes, Vec x, Mat *A, 
                                         Mat *B, MatStructure* flag,
                                         void* ctx)
#else
PetscErrorCode buckettools::FormJacobian(SNES snes, Vec x, Mat A, 
                                         Mat B, void* ctx)
#endif
{
  log(INFO, "In FormJacobian");

  PetscErrorCode perr;                                               // petsc error code

  SNESCtx *snesctx = (SNESCtx *)ctx;                                 // cast the snes context

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  PetscInt iter;
  perr = SNESGetIterationNumber(snes, &iter); CHKERRQ(perr);
  (*solver).iteration_count(iter);

  if ((*solver).monitor_norms())
  {
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormJacobian(1): 2-norm x = %g", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormJacobian(1): inf-norm x = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(1): Frobenius norm A = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*A,NORM_INFINITY,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(1): inf-norm A = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(1): Frobenius norm B = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*B,NORM_INFINITY,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(1): inf-norm B = %g", norm);
  }

  Function_ptr iteratedfunction = (*system).iteratedfunction();      // collect the iterated system bucket function
  const std::vector< const dolfin::DirichletBC* >& bcs = 
                                         (*system).bcs();   // get the vector of bcs
  dolfin::PETScVector iteratedvec(x);
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
  dolfin::PETScMatrix matrix(*A), matrixpc(*B);
  #else
  dolfin::PETScMatrix matrix(A), matrixpc(B);
  #endif

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::SystemAssembler assembler((*solver).bilinear_form(), (*solver).linear_form(),
                                    bcs);
  assembler.assemble(matrix);                                        // assemble the matrix from the context bilinear form
  if ((*solver).ident_zeros())
  {
    matrix.ident_zeros();
  }

  if ((*solver).bilinearpc_form())                                   // do we have a different bilinear pc form associated?
  {
    dolfin::SystemAssembler assemblerpc((*solver).bilinearpc_form(), (*solver).linear_form(),
                                      bcs);
    assemblerpc.assemble(matrixpc);
    if ((*solver).ident_zeros_pc())
    {
      matrixpc.ident_zeros();
    }
  }

  for (Form_const_it f_it = (*solver).solverforms_begin();           // update any solver forms/matrices/submatrices as well
                     f_it != (*solver).solverforms_end(); f_it++)    // - these will already be attached to the appropriate
  {                                                                  // ksps so be careful just to update their pointers
    PETScMatrix_ptr solvermatrix = (*solver).fetch_solvermatrix((*f_it).first);
    dolfin::SystemAssembler assemblerform((*f_it).second, (*solver).linear_form(),
                                      bcs);
    assemblerform.assemble(*solvermatrix);
    if((*solver).solverident_zeros((*f_it).first))
    {
      (*solvermatrix).ident_zeros();
    }

    IS is = (*solver).fetch_solverindexset((*f_it).first);
    Mat submatrix = (*solver).fetch_solversubmatrix((*f_it).first);
    perr = MatGetSubMatrix((*solvermatrix).mat(), is, is, MAT_REUSE_MATRIX, &submatrix);
    CHKERRQ(perr);

  }

  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
  *flag = SAME_NONZERO_PATTERN;                                      // both matrices are assumed to have the same sparsity
  #endif

  if ((*solver).monitor_norms())
  {
    PetscReal norm;

    perr = VecNorm(x,NORM_2,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormJacobian(2): 2-norm x = %g", norm);

    perr = VecNorm(x,NORM_INFINITY,&norm); CHKERRQ(perr);
    log(dolfin::get_log_level(), "FormJacobian(2): inf-norm x = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(2): Frobenius norm A = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*A,NORM_INFINITY,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(2): inf-norm A = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(2): Frobenius norm B = %g", norm);

    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    perr = MatNorm(*B,NORM_INFINITY,&norm); CHKERRQ(perr);
    #else
    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRQ(perr);
    #endif
    log(dolfin::get_log_level(), "FormJacobian(2): inf-norm B = %g", norm);
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc snes monitor callback function that outputs a visualization file
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::SNESCustomMonitor(SNES snes, PetscInt its,
                                      PetscReal norm, void* mctx)
{
  log(INFO, "In SNESCustomMonitor");

  std::stringstream buffer;                                          // string buffer
  PetscErrorCode perr;                                               // petsc error code

  CustomMonitorCtx *snesctx = (CustomMonitorCtx *)mctx;

  SolverBucket* solver = (*snesctx).solver;                          // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  (*solver).iteration_count(its);                                    // set the iteration count

  Vec x;
  perr = SNESGetSolution(snes, &x);  CHKERRQ(perr);                  // get the solution vector from snes
  dolfin::PETScVector sol(x);

  *(*(*system).iteratedfunction()).vector() = sol;

  Vec rx;
  perr = SNESGetFunction(snes, &rx, PETSC_NULL, PETSC_NULL);         // get the residual function vector from snes
  CHKERRQ(perr);
  dolfin::PETScVector solresid(rx);

  *(*(*system).residualfunction()).vector() = solresid;

  Vec dx;
  perr = SNESGetSolutionUpdate(snes, &dx);  CHKERRQ(perr);           // get the solution update vector from snes
  dolfin::PETScVector solupdate(dx);

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
  log(INFO, "In KSPCustomMonitor");

  std::stringstream buffer;                                          // string buffer
  PetscErrorCode perr;                                               // petsc error code

  CustomMonitorCtx *kspctx = (CustomMonitorCtx *)mctx;

  SolverBucket* solver = (*kspctx).solver;                           // retrieve a (standard) pointer to this solver
  SystemBucket* system = (*solver).system();                         // retrieve a (standard) pointer to the parent system of this solver
  Bucket*       bucket = (*system).bucket();                         // retrieve a (standard) pointer to the parent bucket of this solver

  Vec x;
  perr = KSPBuildSolution(ksp, PETSC_NULL, &x);  CHKERRQ(perr);      // get the solution vector from the ksp
  dolfin::PETScVector sol(x);

  *(*(*system).iteratedfunction()).vector() = sol;

  Vec rx;
  perr = KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &rx);  CHKERRQ(perr);// get the residual vector from the ksp
  dolfin::PETScVector solresid(rx);

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
    log(INFO, "In KSPNullSpaceMonitor");

    PetscErrorCode perr;                                             // petsc error code

    MatNullSpace sp;

    perr = KSPGetNullSpace(ksp, &sp); CHKERRQ(perr);

    if (sp)
    {
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
      PetscBool isNull;
      #else
      PetscTruth isNull;
      #endif

      Mat Amat, Pmat;
     
      #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
      MatStructure flag;

      perr = KSPGetOperators(ksp, &Amat, &Pmat, &flag); CHKERRQ(perr);
      #else
      perr = KSPGetOperators(ksp, &Amat, &Pmat); CHKERRQ(perr);
      #endif

      perr = MatNullSpaceTest(sp, Amat, &isNull); CHKERRQ(perr);

      if (!isNull)
      {
        log(WARNING, "MatNullSpaceTest does not believe provided null space is a null space of the matrix.");
      }
      else
      {
        log(INFO, "MatNullSpaceTest thinks provided null space is a null space of the matrix.");
      }
    }

  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// a dummy routine to keep snes vi happy
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu)
{
                                                                     // do nothing
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// check if petsc has failed and throw a sigint if it has
//*******************************************************************|************************************************************//
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
void buckettools::petsc_failure(PetscErrorCode perr,
                                const std::string &filename,
                                const int &line,
                                const std::string &dirname,
                                const std::string &petsc_function)
#else
void buckettools::petsc_failure(PetscErrorCode perr,
                                const std::string &filename,
                                const int &line,
                                const std::string &petsc_function)
#endif
{
  if (PetscUnlikely(perr))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    PetscErrorCode perr2 = PetscError(PETSC_COMM_SELF,line,petsc_function.c_str(),filename.c_str(),dirname.c_str(),perr,PETSC_ERROR_REPEAT," ");
    #else
    PetscErrorCode perr2 = PetscError(PETSC_COMM_SELF,line,petsc_function.c_str(),filename.c_str(),perr,PETSC_ERROR_REPEAT," ");
    #endif

    failure(filename, line, 
            "Call to PETSc function returned an error.",
            "PETSc error code: %d.", perr);
  }
}

//*******************************************************************|************************************************************//
// check if petsc has an error and terminate if it has
//*******************************************************************|************************************************************//
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
void buckettools::petsc_error(PetscErrorCode perr,
                              const std::string &filename,
                              const int &line,
                              const std::string &dirname,
                              const std::string &petsc_function)
#else
void buckettools::petsc_error(PetscErrorCode perr,
                              const std::string &filename,
                              const int &line,
                              const std::string &petsc_function)
#endif
{
  if (PetscUnlikely(perr))
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
    PetscErrorCode perr2 = PetscError(PETSC_COMM_SELF,line,petsc_function.c_str(),filename.c_str(),dirname.c_str(),perr,PETSC_ERROR_REPEAT," ");
    #else
    PetscErrorCode perr2 = PetscError(PETSC_COMM_SELF,line,petsc_function.c_str(),filename.c_str(),perr,PETSC_ERROR_REPEAT," ");
    #endif

    error(filename, line,
          "Call to PETSc function returned an error.",
          "PETSc error code: %d.", perr);
  }
}


