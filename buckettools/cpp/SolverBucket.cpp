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


#include "BoostTypes.h"
#include "SolverBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include "Logger.h"
#include <dolfin.h>
#include <string>
#include <signal.h>

using namespace buckettools;

//*******************************************************************|************************************************************//
// default constructor
//*******************************************************************|************************************************************//
SolverBucket::SolverBucket()
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SolverBucket::SolverBucket(SystemBucket* system) : system_(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SolverBucket::~SolverBucket()
{
  PetscErrorCode perr;                                               // petsc error code

  if(type()=="SNES")
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = SNESDestroy(&snes_); petsc_err(perr);                       // destroy the snes object
    #else
    perr = SNESDestroy(snes_); petsc_err(perr);                        // destroy the snes object
    #endif
  }

  if(type()=="Picard")
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = KSPDestroy(&ksp_); petsc_err(perr);                         // destroy the ksp object
    #else
    perr = KSPDestroy(ksp_); petsc_err(perr);                          // destroy the ksp object
    #endif
  }

  if(sp_)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatNullSpaceDestroy(&sp_); petsc_err(perr);                 // destroy the null space object
    #else
    perr = MatNullSpaceDestroy(sp_); petsc_err(perr);                  // destroy the null space object
    #endif
  }

  for (std::map<std::string, IS>::iterator is_it = solverindexsets_.begin();
                                           is_it != solverindexsets_.end();
                                           is_it++)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = ISDestroy(&(*is_it).second); petsc_err(perr);                  // destroy the IS, necessary?
    #else
    perr = ISDestroy((*is_it).second); petsc_err(perr);                   // destroy the IS, necessary?
    #endif
  }

  for (std::map<std::string, Mat>::iterator mat_it = solversubmatrices_.begin();
                                            mat_it != solversubmatrices_.end();
                                            mat_it++)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = MatDestroy(&(*mat_it).second); petsc_err(perr);
    #else
    perr = MatDestroy((*mat_it).second); petsc_err(perr);
    #endif
  }

  for (ConvergenceFile_it f_it = convergencefiles_begin(); f_it != convergencefiles_end(); f_it++)
  {
    (*(*f_it).second).close();
  }

  for (KSPConvergenceFile_it k_it = kspconvergencefiles_begin(); k_it != kspconvergencefiles_end(); k_it++)
  {
    (*(*k_it).second).close();
  }

}

//*******************************************************************|************************************************************//
// solve the bilinear system described by the forms in the solver bucket
//*******************************************************************|************************************************************//
bool SolverBucket::solve()
{
  PetscErrorCode perr;

  if (solve_location()==SOLVE_NEVER)
  {
    tf_err("Unable to solve as solve_location is set to never.", "SolverBucket name: %s, SystemBucket name: %s",
           name_.c_str(), (*system_).name().c_str());
  }

  log(INFO, "Solving for %s::%s using %s", 
                          (*system_).name().c_str(), name().c_str(), 
                          type().c_str());

  *iteration_count_ = 0;                                             // an iteration counter

  if (type()=="SNES")                                                // this is a petsc snes solver - FIXME: switch to an enumerated type
  {                                                                  // loop over the collected vector of system bcs
    for(std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator    
                      bc = (*system_).bcs_begin(); 
                      bc != (*system_).bcs_end(); bc++)
    {
      (*(*bc)).apply(*(*(*system_).function()).vector());            // apply the bcs to the solution and
      (*(*bc)).apply(*(*(*system_).iteratedfunction()).vector());    // iterated solution
    }
    *work_ = (*(*(*system_).function()).vector());                   // set the work vector to the function vector
    perr = SNESSolve(snes_, PETSC_NULL, (*work_).vec());             // call petsc to perform a snes solve
    petsc_fail(perr);
    snes_check_convergence_();
    (*(*(*system_).function()).vector()) = *work_;                   // update the function
  }
  else if (type()=="Picard")                                         // this is a hand-rolled picard iteration - FIXME: switch to enum
  {


    assert(residual_);                                               // we need to assemble the residual again here as it may depend
                                                                     // on other systems that have been solved since the last call
    dolfin::Assembler assemblerres;
    assemblerres.assemble(*res_, *residual_);                        // assemble the residual
    for(std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator bc = 
                          (*system_).bcs_begin(); 
                          bc != (*system_).bcs_end(); bc++)
    {                                                                // apply bcs to residuall (should we do this?!)
      (*(*bc)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
    }

    double aerror = (*res_).norm("l2");                              // work out the initial absolute l2 error (this should be
                                                                     // initialized to the right value on the first pass and still
                                                                     // be the correct value from the previous sweep (if this stops
                                                                     // being the case it will be necessary to assemble the residual
                                                                     // here too)
    double aerror0 = aerror;                                         // record the initial absolute error
    double rerror;
    if(aerror==0.0)
    {
      rerror = aerror;                                               // relative error, starts out as 0 - won't iterate at all
    }
    else
    {
      rerror = aerror/aerror0;                                       // relative error, starts out as 1.
    }

    log(INFO, "  %u Picard Residual Norm (absolute, relative) = %g, %g\n", 
                                    iteration_count(), aerror, rerror);

    ConvergenceFile_ptr convfile = convergencefile();
    if(convfile)
    {
      *(*(*system()).residualfunction()).vector() = (*std::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));
      if (convfile)
      {
        (*convfile).write_data();
      }
    }


    (*(*(*system_).iteratedfunction()).vector()) =                   // system iterated function gets set to the function values
                                (*(*(*system_).function()).vector());

    while (iteration_count() < minits_ ||                            // loop for the minimum number of iterations or
          (iteration_count() < maxits_ &&                            // up to the maximum number of iterations 
                           rerror > rtol_ && aerror > atol_))        // until the max is reached or a tolerance criterion is
    {                                                                // satisfied
      (*iteration_count_)++;                                         // increment iteration counter

      dolfin::SystemAssembler assembler(bilinear_, linear_,
                                        (*system_).bcs());
      assembler.assemble(*matrix_, *rhs_);

      if(ident_zeros_)
      {
        (*matrix_).ident_zeros();
      }

      if (bilinearpc_)                                               // if there's a pc associated
      {
        assert(matrixpc_);
        dolfin::SystemAssembler assemblerpc(bilinearpc_, linear_,
                                          (*system_).bcs());
        assemblerpc.assemble(*matrixpc_);

        if(ident_zeros_pc_)
        {
          (*matrixpc_).ident_zeros();
        }

        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
        perr = KSPSetOperators(ksp_, (*matrix_).mat(),              // set the ksp operators with two matrices
                                     (*matrixpc_).mat(), 
                                     SAME_NONZERO_PATTERN); 
        #else
        perr = KSPSetOperators(ksp_, (*matrix_).mat(),              // set the ksp operators with two matrices
                                     (*matrixpc_).mat()); 
        #endif
        petsc_err(perr);
      }
      else
      {
        #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 5
        perr = KSPSetOperators(ksp_, (*matrix_).mat(),              // set the ksp operators with the same matrices
                                      (*matrix_).mat(), 
                                        SAME_NONZERO_PATTERN); 
        #else
        perr = KSPSetOperators(ksp_, (*matrix_).mat(),              // set the ksp operators with the same matrices
                                      (*matrix_).mat()); 
        #endif
        petsc_err(perr);
      }

      for (Form_const_it f_it = solverforms_begin(); 
                         f_it != solverforms_end(); f_it++)
      {
        PETScMatrix_ptr solvermatrix = solvermatrices_[(*f_it).first];
        dolfin::SystemAssembler assemblerform((*f_it).second, linear_,
                                          (*system_).bcs());
        assemblerform.assemble(*solvermatrix);

        if(solverident_zeros_[(*f_it).first])
        {
          (*solvermatrix).ident_zeros();
        }

        IS is = solverindexsets_[(*f_it).first];
        Mat submatrix = solversubmatrices_[(*f_it).first];
        perr = MatGetSubMatrix((*solvermatrix).mat(), is, is, MAT_REUSE_MATRIX, &submatrix);
        petsc_err(perr);

      }

      if (monitor_norms())
      {
        PetscReal norm;

        perr = VecNorm((*rhs_).vec(),NORM_2,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: 2-norm rhs = %g", norm);

        perr = VecNorm((*rhs_).vec(),NORM_INFINITY,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: inf-norm rhs = %g", norm);

        perr = VecNorm((*work_).vec(),NORM_2,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: 2-norm work = %g", norm);

        perr = VecNorm((*work_).vec(),NORM_INFINITY,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: inf-norm work = %g", norm);

        perr = MatNorm((*matrix_).mat(),NORM_FROBENIUS,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: Frobenius norm matrix = %g", norm);

        perr = MatNorm((*matrix_).mat(),NORM_INFINITY,&norm); petsc_err(perr);
        log(dolfin::get_log_level(), "Picard: inf-norm matrix = %g", norm);

        if (bilinearpc_)
        {
          perr = MatNorm((*matrixpc_).mat(),NORM_FROBENIUS,&norm); petsc_err(perr);
          log(dolfin::get_log_level(), "Picard: Frobenius norm matrix pc = %g", norm);

          perr = MatNorm((*matrixpc_).mat(),NORM_INFINITY,&norm); petsc_err(perr);
          log(dolfin::get_log_level(), "Picard: inf-norm matrix pc = %g", norm);
        }
      }

      *work_ = (*(*(*system_).iteratedfunction()).vector());         // set the work vector to the iterated function
      perr = KSPSolve(ksp_, (*rhs_).vec(), (*work_).vec());        // perform a linear solve
      petsc_fail(perr);
      ksp_check_convergence_(ksp_);
      (*(*(*system_).iteratedfunction()).vector()) = *work_;         // update the iterated function with the work vector

      assert(residual_);
      assemblerres.assemble(*res_, *residual_);                      // assemble the residual
      for(std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator bc = 
                             (*system_).bcs_begin(); 
                             bc != (*system_).bcs_end(); bc++)
      {                                                              // apply bcs to residual (should we do this?!)
        (*(*bc)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
      }

      aerror = (*res_).norm("l2");                                   // work out absolute error
      rerror = aerror/aerror0;                                       // and relative error
      log(INFO, "  %u Picard Residual Norm (absolute, relative) = %g, %g\n", 
                          iteration_count(), aerror, rerror);
                                                                     // and decide to loop or not...

      if(convfile)
      {
        *(*(*system()).residualfunction()).vector() = (*std::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));
        if (convfile)
        {
          (*convfile).write_data();
        }
      }

    }

    if (iteration_count() == maxits_ && rerror > rtol_ && aerror > atol_)
    {
      log(WARNING, "it = %d, maxits_ = %d", iteration_count(), maxits_);
      log(WARNING, "rerror = %.12e, rtol_ = %.12e", rerror, rtol_);
      log(WARNING, "aerror = %.12e, atol_ = %.12e", aerror, atol_);
      if (ignore_failures_)
      {
        log(WARNING, "Ignoring: Picard failure. Solver: %s::%s.", (*system_).name().c_str(), name().c_str());
      }
      else
      {
        tf_fail("Picard iterations failed to converge.", "Iteration count, relative error or absolute error too high.");
      }
    }

    (*(*(*system_).function()).vector()) =                              // update the function values with the iterated values
                      (*(*(*system_).iteratedfunction()).vector());

  }
  else                                                               // don't know what solver type this is
  {
    tf_err("Unknown solver type.", "Type: %s", type_.c_str());
  }

  (*system()).postprocess_values();

  *(*(*system()).residualfunction()).vector() = (*std::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));

  if (solved_)
  {
    *solved_ = true;
  }

  return true;

}

//*******************************************************************|************************************************************//
// return the l2 norm of the residual
//*******************************************************************|************************************************************//
double SolverBucket::residual_norm()
{
  assert(residual_);
  dolfin::Assembler assembler;

  assembler.assemble(*res_, *residual_);
  for(std::vector< std::shared_ptr<const dolfin::DirichletBC> >::const_iterator bc = 
                        (*system_).bcs_begin(); 
                        bc != (*system_).bcs_end(); bc++)
  {                                                                  // apply bcs to residual (should we do this?!)
    (*(*bc)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
  }

  *(*(*system()).residualfunction()).vector() = (*std::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));

  double norm  = (*res_).norm("l2");
  return norm;
}

//*******************************************************************|************************************************************//
// update the solver at the end of a timestep
//*******************************************************************|************************************************************//
void SolverBucket::resetcalculated()
{
  if (solved_)
  {
    *solved_ = false;
  }
  *iteration_count_ = 0;                                             // an iteration counter
}

//*******************************************************************|************************************************************//
// loop over the forms in this solver bucket and attach the coefficients they request using the parent bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::attach_form_coeffs()
{
  (*(*system_).bucket()).attach_coeffs(forms_begin(), forms_end());
}

//*******************************************************************|************************************************************//
// create a null space object
//*******************************************************************|************************************************************//
void SolverBucket::create_nullspace()
{
  std::size_t nnulls = nullspacevectors_.size();
  if (nnulls > 0)
  {
    PetscErrorCode perr;                                             // petsc error code
    Vec vecs[nnulls];
    for (uint i = 0; i < nnulls; i++)
    {
      vecs[i] = (*(nullspacevectors_[i])).vec();
    }
    orthonormalize_petsc_vecs_(vecs, nnulls);
    perr = MatNullSpaceCreate((*(nullspacevectors_[0])).mpi_comm(), 
                              PETSC_FALSE, nnulls, vecs, &sp_); 
    petsc_err(perr);
  }
  else
  {
    sp_ = PETSC_NULL;
  }
}

//*******************************************************************|************************************************************//
// return the number of nonlinear iterations taken
//*******************************************************************|************************************************************//
const int SolverBucket::iteration_count() const
{
  return *iteration_count_;
}

//*******************************************************************|************************************************************//
// set the number of nonlinear iterations taken
//*******************************************************************|************************************************************//
void SolverBucket::iteration_count(const int &it)
{
  *iteration_count_ = it;
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_form(Form_ptr form, const std::string &name)
{
  Form_hash_it f_it = forms_.get<om_key_hash>().find(name);                                  // check if this name already exists
  if (f_it != forms_.get<om_key_hash>().end())
  {
    tf_err("Form already exists in solver.", "Form name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    forms_.insert(om_item<const std::string, Form_ptr>(name, form));                                             // if not, register the form in the maps
  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the solver bucket contains a form with the given name
//*******************************************************************|************************************************************//
bool SolverBucket::contains_form(const std::string &name)                   
{
  Form_hash_it f_it = forms_.get<om_key_hash>().find(name);
  return f_it != forms_.get<om_key_hash>().end();
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a form from the solver bucket data maps
//*******************************************************************|************************************************************//
Form_ptr SolverBucket::fetch_form(const std::string &name)
{
  Form_hash_it f_it = forms_.get<om_key_hash>().find(name);                                  // check if this name already exists
  if (f_it == forms_.get<om_key_hash>().end())
  {
    tf_err("Form does not exist in solver.", "Form name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the forms_ map
//*******************************************************************|************************************************************//
Form_it SolverBucket::forms_begin()
{
  return forms_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the forms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::forms_begin() const
{
  return forms_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the forms_ map
//*******************************************************************|************************************************************//
Form_it SolverBucket::forms_end()
{
  return forms_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the forms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::forms_end() const
{
  return forms_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_solverform(Form_ptr form, const std::string &name)
{
  Form_hash_it f_it = solverforms_.get<om_key_hash>().find(name);                            // check if this name already exists
  if (f_it != solverforms_.get<om_key_hash>().end())
  {
    tf_err("Solver form already exists in solver.", "Solver form name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    solverforms_.insert(om_item<const std::string, Form_ptr>(name, form));                                             // if not, register the form in the maps
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the solverforms_ map
//*******************************************************************|************************************************************//
Form_it SolverBucket::solverforms_begin()
{
  return solverforms_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the solverforms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::solverforms_begin() const
{
  return solverforms_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the solverforms_ map
//*******************************************************************|************************************************************//
Form_it SolverBucket::solverforms_end()
{
  return solverforms_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the solverforms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::solverforms_end() const
{
  return solverforms_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a petsc matrix from the solver bucket data maps
//*******************************************************************|************************************************************//
PETScMatrix_ptr SolverBucket::fetch_solvermatrix(const std::string &name)
{
  std::map< std::string, PETScMatrix_ptr>::iterator s_it = 
                                          solvermatrices_.find(name);// check if this name already exists
  if (s_it == solvermatrices_.end())
  {
    tf_err("Solver matrix does not exist in solver.", "Solver matrix name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to an index set for this solver sub matrix
//*******************************************************************|************************************************************//
IS SolverBucket::fetch_solverindexset(const std::string &name)
{
  std::map< std::string, IS >::iterator s_it = 
                                         solverindexsets_.find(name);// check if this name already exists
  if (s_it == solverindexsets_.end())
  {
    tf_err("Solver index set does not exist in solver.", "Solver index set name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a bool indicating if the named solver form/matrix should have zeros idented
//*******************************************************************|************************************************************//
bool SolverBucket::solverident_zeros(const std::string &name)
{
  std::map< std::string, bool >::iterator s_it = 
                                       solverident_zeros_.find(name);// check if this name already exists
  if (s_it == solverident_zeros_.end())
  {
    tf_err("Solver ident zeros does not exist in solver.", "Solver ident zeros name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a petsc matrix from the solver bucket data maps
//*******************************************************************|************************************************************//
Mat SolverBucket::fetch_solversubmatrix(const std::string &name)
{
  std::map< std::string, Mat >::iterator s_it = 
                                          solversubmatrices_.find(name);// check if this name already exists
  if (s_it == solversubmatrices_.end())
  {
    tf_err("Solver sub matrix does not exist in solver.", "Solver sub matrix name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// register a pointer to a systems solver in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_systemssolver(SystemsSolverBucket* solver, const std::string &name)
{
  p_SystemsSolverBucket_hash_it s_it = systemssolvers_.get<om_key_hash>().find(name);
  if (s_it == systemssolvers_.end())
  {
    systemssolvers_.insert(om_item<const std::string, SystemsSolverBucket*>(name, solver));
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a systemssolver from the solver bucket data maps
//*******************************************************************|************************************************************//
SystemsSolverBucket* SolverBucket::fetch_systemssolver(const std::string &name)
{
  p_SystemsSolverBucket_hash_it s_it = systemssolvers_.get<om_key_hash>().find(name);                                  // check if this name already exists
  if (s_it == systemssolvers_.get<om_key_hash>().end())
  {
    tf_err("SystemsSolverBucket does not exist in solver.", "SystemsSolverBucket name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*s_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the systemssolvers_ map
//*******************************************************************|************************************************************//
p_SystemsSolverBucket_it SolverBucket::systemssolvers_begin()
{
  return systemssolvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the systemssolvers_ map
//*******************************************************************|************************************************************//
p_SystemsSolverBucket_const_it SolverBucket::systemssolvers_begin() const
{
  return systemssolvers_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the systemssolvers_ map
//*******************************************************************|************************************************************//
p_SystemsSolverBucket_it SolverBucket::systemssolvers_end()
{
  return systemssolvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the systemssolvers_ map
//*******************************************************************|************************************************************//
p_SystemsSolverBucket_const_it SolverBucket::systemssolvers_end() const
{
  return systemssolvers_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a pointer to a convergence file in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_convergencefile(ConvergenceFile_ptr convfile, const std::string &name)
{
  ConvergenceFile_hash_it f_it = convergencefiles_.get<om_key_hash>().find(name);
  if (f_it != convergencefiles_.get<om_key_hash>().end())
  {
    tf_err("ConvergenceFile already exists in solver.", "ConvergenceFile name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    convergencefiles_.insert(om_item<const std::string, ConvergenceFile_ptr>(name, convfile));
  }
}

//*******************************************************************|************************************************************//
// return a pointer to the convergencefile for the current systems solver (if one exists)
//*******************************************************************|************************************************************//
ConvergenceFile_ptr SolverBucket::convergencefile()
{
  ConvergenceFile_hash_it f_it = convergencefiles_.get<om_key_hash>().find(current_systemssolver());  // check if this name already exists
  if (f_it == convergencefiles_.get<om_key_hash>().end())
  {
    return NULL;
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a convergencefile from the solver bucket data maps
//*******************************************************************|************************************************************//
ConvergenceFile_ptr SolverBucket::fetch_convergencefile(const std::string &name)
{
  ConvergenceFile_hash_it f_it = convergencefiles_.get<om_key_hash>().find(name);                                  // check if this name already exists
  if (f_it == convergencefiles_.get<om_key_hash>().end())
  {
    tf_err("ConvergenceFile does not exist in solver.", "ConvergenceFile name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the convergencefiles_ map
//*******************************************************************|************************************************************//
ConvergenceFile_it SolverBucket::convergencefiles_begin()
{
  return convergencefiles_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the convergencefiles_ map
//*******************************************************************|************************************************************//
ConvergenceFile_const_it SolverBucket::convergencefiles_begin() const
{
  return convergencefiles_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the convergencefiles_ map
//*******************************************************************|************************************************************//
ConvergenceFile_it SolverBucket::convergencefiles_end()
{
  return convergencefiles_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the convergencefiles_ map
//*******************************************************************|************************************************************//
ConvergenceFile_const_it SolverBucket::convergencefiles_end() const
{
  return convergencefiles_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// register a pointer to a ksp convergence file in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_kspconvergencefile(KSPConvergenceFile_ptr kspconvfile, const std::string &name)
{
  KSPConvergenceFile_hash_it f_it = kspconvergencefiles_.get<om_key_hash>().find(name);
  if (f_it != kspconvergencefiles_.get<om_key_hash>().end())
  {
    tf_err("ConvergenceFile already exists in solver.", "ConvergenceFile name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    kspconvergencefiles_.insert(om_item<const std::string, KSPConvergenceFile_ptr>(name, kspconvfile));
  }
}

//*******************************************************************|************************************************************//
// return a pointer to a ksp convergencefile from the solver bucket data maps
//*******************************************************************|************************************************************//
KSPConvergenceFile_ptr SolverBucket::fetch_kspconvergencefile(const std::string &name)
{
  KSPConvergenceFile_hash_it f_it = kspconvergencefiles_.get<om_key_hash>().find(name);                                  // check if this name already exists
  if (f_it == kspconvergencefiles_.get<om_key_hash>().end())
  {
    tf_err("KSPConvergenceFile does not exist in solver.", "KSPConvergenceFile name: %s, SolverBucket name: %s, SystemBucket name: %s",
           name.c_str(), name_.c_str(), (*system_).name().c_str());
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return a pointer to the kspconvergencefile for the current systems solver (if one exists)
//*******************************************************************|************************************************************//
KSPConvergenceFile_ptr SolverBucket::kspconvergencefile()
{
  KSPConvergenceFile_hash_it f_it = kspconvergencefiles_.get<om_key_hash>().find(current_systemssolver());  // check if this name already exists
  if (f_it == kspconvergencefiles_.get<om_key_hash>().end())
  {
    return NULL;
  }
  else
  {
    return (*f_it).second;                                           // if it does, return it
  }
}

//*******************************************************************|************************************************************//
// return an iterator to the beginning of the kspconvergencefiles_ map
//*******************************************************************|************************************************************//
KSPConvergenceFile_it SolverBucket::kspconvergencefiles_begin()
{
  return kspconvergencefiles_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the kspconvergencefiles_ map
//*******************************************************************|************************************************************//
KSPConvergenceFile_const_it SolverBucket::kspconvergencefiles_begin() const
{
  return kspconvergencefiles_.get<om_key_seq>().begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the kspconvergencefiles_ map
//*******************************************************************|************************************************************//
KSPConvergenceFile_it SolverBucket::kspconvergencefiles_end()
{
  return kspconvergencefiles_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the kspconvergencefiles_ map
//*******************************************************************|************************************************************//
KSPConvergenceFile_const_it SolverBucket::kspconvergencefiles_end() const
{
  return kspconvergencefiles_.get<om_key_seq>().end();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the solver bucket
//*******************************************************************|************************************************************//
const std::string SolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << std::endl;
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the forms in the solver bucket
//*******************************************************************|************************************************************//
const std::string SolverBucket::forms_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = forms_begin(); f_it != forms_end(); f_it++ )
  {
    s << indentation << "Form " << (*f_it).first  << std::endl;
  }
  return s.str();
}

//*******************************************************************|************************************************************//
// report the convergence of the snes solver
//*******************************************************************|************************************************************//
void SolverBucket::snes_check_convergence_()
{
  PetscErrorCode perr;                                               // petsc error code

  assert(type()=="SNES");

  log(INFO, "Convergence for %s::%s", 
                          (*system_).name().c_str(), name().c_str());

  SNESConvergedReason snesreason;                                    // check what the convergence reason was
  PetscInt snesiterations;
  PetscInt sneslsiterations;
//  const char **snesprefix;
//  perr = SNESGetOptionsPrefix(snes_, snesprefix); petsc_err(perr);   // FIXME: segfaults!
  perr = SNESGetConvergedReason(snes_, &snesreason); petsc_err(perr);     
  perr = SNESGetIterationNumber(snes_, &snesiterations); petsc_err(perr);
  iteration_count(snesiterations);
  perr = SNESGetLinearSolveIterations(snes_, &sneslsiterations);  petsc_err(perr);
  log(INFO, "SNESConvergedReason %d", snesreason);
  log(INFO, "SNES n/o iterations %d", 
                              snesiterations);
  log(INFO, "SNES n/o linear solver iterations %d", 
                              sneslsiterations);
  if (snesreason<0)
  {
    if (ignore_failures_)
    {
      log(WARNING, "Ignoring: SNES failure. Solver: %s::%s, SNESConvergedReason %d.", (*system_).name().c_str(), name().c_str(), snesreason);
    }
    else
    {
      tf_fail("SNES failed to converge.", "Solver: %s::%s, SNESConvergedReason: %d.", (*system_).name().c_str(), name().c_str(), snesreason);
    }
  }

  ksp_check_convergence_(ksp_, 1);

}

//*******************************************************************|************************************************************//
// report the convergence of a ksp solver
//*******************************************************************|************************************************************//
void SolverBucket::ksp_check_convergence_(KSP &ksp, int indent)
{
  PetscErrorCode perr;                                               // petsc error code
  std::string indentation (indent*2, ' ');

  if (indent==0)
  {
    log(INFO, "Convergence for %s::%s", 
                          (*system_).name().c_str(), name().c_str());
  }

  KSPConvergedReason kspreason;                                      // check what the convergence reason was
  PetscInt kspiterations;
  //const char **kspprefix;
  //perr = KSPGetOptionsPrefix(ksp, kspprefix); petsc_err(perr);       // FIXME: segfaults!
  perr = KSPGetConvergedReason(ksp, &kspreason); petsc_err(perr);     
  perr = KSPGetIterationNumber(ksp, &kspiterations); petsc_err(perr);     
  log(INFO, "%sKSPConvergedReason %d", 
                              indentation.c_str(), kspreason);
  log(INFO, "%sKSP n/o iterations %d", 
                              indentation.c_str(), kspiterations);
  if (indent==0 && kspreason<0)
  {
    if (ignore_failures_)
    {
      log(WARNING, "Ignoring: KSP failure. Solver: %s::%s, KSPConvergedReason %d.", (*system_).name().c_str(), name().c_str(), kspreason);
    }
    else
    {
      tf_fail("KSP failed to converge.", "Solver: %s::%s, KSPConvergedReason: %d.", (*system_).name().c_str(), name().c_str(), kspreason);
    }
  }


  indent++;

  PC pc;
  perr = KSPGetPC(ksp, &pc); petsc_err(perr);
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 3
  PCType pctype;
  #else
  const PCType pctype;
  #endif
  perr = PCGetType(pc, &pctype); petsc_err(perr);

  if ((std::string)pctype=="ksp")
  {
    KSP subksp;                                                      // get the subksp from this pc
    perr = PCKSPGetKSP(pc, &subksp); petsc_err(perr);
    ksp_check_convergence_(subksp, indent);
  }
  else if ((std::string)pctype=="fieldsplit")
  {
    KSP *subksps;                                                    // get the fieldsplit subksps
    PetscInt nsubksps;
    perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); 
    petsc_err(perr); 
    for (PetscInt i = 0; i < nsubksps; i++)
    {
      ksp_check_convergence_(subksps[i], indent);
    }
  }
  
}

//*******************************************************************|************************************************************//
// checkpoint the solverbucket
//*******************************************************************|************************************************************//
void SolverBucket::checkpoint()
{
  checkpoint_options_();
}

//*******************************************************************|************************************************************//
// orthonormalize an array of petsc Vecs
//*******************************************************************|************************************************************//
void SolverBucket::orthonormalize_petsc_vecs_(Vec vecs[], PetscInt n)
{
  PetscErrorCode perr;                                               // petsc error code
  PetscScalar *dots;
  perr = PetscMalloc1(n-1,&dots); petsc_err(perr);
  perr = VecNormalize(vecs[0], PETSC_NULL); petsc_err(perr);
  for (PetscInt i=1; i<n; i++) 
  {
    perr = VecMDot(vecs[i],i,vecs,dots); petsc_err(perr);
    for (PetscInt j=0; j<i; j++) 
    {
      dots[j] = -dots[j];
    }
    perr = VecMAXPY(vecs[i],i,dots,vecs); petsc_err(perr);
    perr = VecNormalize(vecs[i], PETSC_NULL); petsc_err(perr);
  }
}

//*******************************************************************|************************************************************//
// virtual checkpointing of options
//*******************************************************************|************************************************************//
void SolverBucket::checkpoint_options_()
{
  tf_err("Failed to find virtual function checkpoint_options_.", "Need to implement a checkpointing method.");
}


