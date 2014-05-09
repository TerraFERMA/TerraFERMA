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


#include "BucketPETScBase.h"
#include <dolfin.h>
#include "petscsnes.h"
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
  dolfin::PETScVector rhs(f), iteratedvec(x);

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::Assembler assembler;
  assembler.assemble(rhs, *(*solver).linear_form());
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
  dolfin::PETScVector iteratedvec(x);
  dolfin::PETScMatrix matrix(*A), matrixpc(*B);

  (*(*iteratedfunction).vector()) = iteratedvec;                     // update the iterated system bucket function

  (*bucket).update_nonlinear();                                      // update nonlinear coefficients

  dolfin::SystemAssembler assembler((*solver).bilinear_form(), (*solver).linear_form(),
                                    bcs);
  assembler.assemble(matrix);                                        // assemble the matrix from the context bilinear form
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
    dolfin::SystemAssembler assemblerpc((*solver).bilinearpc_form(), (*solver).linear_form(),
                                      bcs);
    assemblerpc.assemble(matrixpc);
    for(uint i = 0; i < points.size(); ++i)                          // loop over the points
    {
      (*points[i]).apply(matrixpc);
    }
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
    for(uint i = 0; i < points.size(); ++i)                          // loop over the points
    {
      (*points[i]).apply(*solvermatrix);
    }
    if((*solver).solverident_zeros((*f_it).first))
    {
      (*solvermatrix).ident_zeros();
    }

    IS_ptr is = (*solver).fetch_solverindexset((*f_it).first);
    Mat_ptr submatrix = (*solver).fetch_solversubmatrix((*f_it).first);
    perr = MatGetSubMatrix((*solvermatrix).mat(), *is, *is, MAT_REUSE_MATRIX, &(*submatrix));
    CHKERRQ(perr);

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
// define the petsc snes monitor callback function that outputs a visualization file and a convergence file
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::SNESCustomMonitor(SNES snes, PetscInt its,
                                      PetscReal norm, void* mctx)
{
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

  std::vector< GenericFunction_ptr > functions;                      // create a list of subfunctions
  for (FunctionBucket_const_it f_it = (*system).fields_begin(); 
                               f_it != (*system).fields_end(); 
                                                        f_it++)
  {
    functions.push_back((*(*f_it).second).iteratedfunction());
    functions.push_back((*(*f_it).second).residualfunction());
    functions.push_back((*(*f_it).second).snesupdatefunction());
  }

  if (((*solver).visualization_monitor()))
  {
    if (its==0)
    {
      buffer.str(""); buffer << (*bucket).output_basename() << "_" 
                             << (*system).name() << "_" 
                             << (*solver).name() << "_" 
                             << (*bucket).timestep_count() << "_" 
                             << (*bucket).iteration_count() << "_snes.pvd";
      (*snesctx).pvdfile.reset( new dolfin::File(buffer.str(), "compressed") );
    }
    assert((*snesctx).pvdfile);

    Mesh_ptr sysmesh = (*system).mesh();
    FunctionSpace_ptr visfuncspace = (*bucket).fetch_visfunctionspace(sysmesh);
    (*(*snesctx).pvdfile).write(functions, *visfuncspace, (double) its);
  }

  if ((*solver).convergence_file())
  {

    (*(*solver).convergence_file()).write_data();

  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ksp monitor callback function that outputs a visualization file
//*******************************************************************|************************************************************//
PetscErrorCode buckettools::KSPCustomMonitor(KSP ksp, int it,
                                      PetscReal rnorm, void* mctx)
{
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

  std::vector< GenericFunction_ptr > functions;                      // create a list of subfunctions
  for (FunctionBucket_const_it f_it = (*system).fields_begin(); 
                               f_it != (*system).fields_end(); 
                                                        f_it++)
  {
    functions.push_back((*(*f_it).second).iteratedfunction());
    functions.push_back((*(*f_it).second).residualfunction());
  }

  if ((*solver).type()=="SNES")
  {
    PetscInt iter;
    perr = SNESGetIterationNumber((*solver).snes(), &iter); CHKERRQ(perr);
    (*solver).iteration_count(iter);
  }

  if (((*solver).kspvisualization_monitor()))
  {
    if (it==0)
    {
      buffer.str(""); buffer << (*bucket).output_basename() << "_" 
                             << (*system).name() << "_" 
                             << (*solver).name() << "_" 
                             << (*bucket).timestep_count() << "_" 
                             << (*bucket).iteration_count() << "_"
                             << (*solver).iteration_count() << "_ksp.pvd";
      (*kspctx).pvdfile.reset( new dolfin::File(buffer.str(), "compressed") );
    }
    assert((*kspctx).pvdfile);

    Mesh_ptr sysmesh = (*system).mesh();
    FunctionSpace_ptr visfuncspace = (*bucket).fetch_visfunctionspace(sysmesh);
    (*(*kspctx).pvdfile).write(functions, *visfuncspace, (double) it);
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


PetscErrorCode buckettools::SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu)
{
                                                                     // do nothing
  PetscFunctionReturn(0);
}


