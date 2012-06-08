
#include "BoostTypes.h"
#include "SolverBucket.h"
#include "SystemBucket.h"
#include "Bucket.h"
#include "SignalHandler.h"
#include <dolfin.h>
#include <string>
#include <signal.h>
#include "BucketDolfinBase.h"

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
  empty_();                                                          // empty the solver bucket data structures

  PetscErrorCode perr;                                               // petsc error code

  if(type()=="SNES" && !copy_)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = SNESDestroy(&snes_); CHKERRV(perr);                        // destroy the snes object
    #else
    perr = SNESDestroy(snes_); CHKERRV(perr);                         // destroy the snes object
    #endif
  }

  if(type()=="Picard" && !copy_)
  {
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR > 1
    perr = KSPDestroy(&ksp_); CHKERRV(perr);                          // destroy the ksp object
    #else
    perr = KSPDestroy(ksp_); CHKERRV(perr);                          // destroy the ksp object
    #endif
  }

}

//*******************************************************************|************************************************************//
// solve the bilinear system described by the forms in the solver bucket
//*******************************************************************|************************************************************//
void SolverBucket::solve()
{
  PetscErrorCode perr;

  dolfin::log(dolfin::INFO, "Solving for %s::%s using %s", 
                          (*system_).name().c_str(), name().c_str(), 
                          type().c_str());

  if (type()=="SNES")                                                // this is a petsc snes solver - FIXME: switch to an enumerated type
  {
    for(std::vector< const dolfin::DirichletBC* >::const_iterator    // loop over the collected vector of system bcs
                      bc = (*system_).dirichletbcs_begin(); 
                      bc != (*system_).dirichletbcs_end(); bc++)
    {
      (*(*bc)).apply(*(*(*system_).function()).vector());            // apply the bcs to the solution and
      (*(*bc)).apply(*(*(*system_).iteratedfunction()).vector());    // iterated solution
    }
    *work_ = (*(*(*system_).function()).vector());                   // set the work vector to the function vector
    perr = SNESSolve(snes_, PETSC_NULL, *(*work_).vec());            // call petsc to perform a snes solve
    CHKERRV(perr);
    snes_check_convergence_();
    (*(*(*system_).function()).vector()) = *work_;                   // update the function
  }
  else if (type()=="Picard")                                         // this is a hand-rolled picard iteration - FIXME: switch to enum
  {

    *iteration_count_ = 0;                                           // an iteration counter

    assert(residual_);                                               // we need to assemble the residual again here as it may depend
                                                                     // on other systems that have been solved since the last call
    dolfin::assemble(*res_, *residual_, false);                      // assemble the residual
    for(std::vector< const dolfin::DirichletBC* >::const_iterator bc = 
                          (*system_).dirichletbcs_begin(); 
                          bc != (*system_).dirichletbcs_end(); bc++)
    {                                                                // apply bcs to residuall (should we do this?!)
      (*(*bc)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
    }

    for(std::vector<ReferencePoints_ptr>::const_iterator p = 
                                    (*system_).points_begin(); 
                                p != (*system_).points_end(); p++)
    {                                                                // apply reference points to residual (should we do this?!)
      (*(*p)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
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

    dolfin::info("  %u Picard Residual Norm (absolute, relative) = %g, %g\n", 
                                    iteration_count(), aerror, rerror);

    if(convfile_)
    {
      *(*(*system()).residualfunction()).vector() = (*boost::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));
      if (convfile_)
      {
        (*convfile_).write_data();
      }
    }


    (*(*(*system_).iteratedfunction()).vector()) =                   // system iterated function gets set to the function values
                                (*(*(*system_).function()).vector());

    while (iteration_count() < minits_ ||                            // loop for the minimum number of iterations or
          (iteration_count() < maxits_ &&                            // up to the maximum number of iterations 
                           rerror > rtol_ && aerror > atol_))        // until the max is reached or a tolerance criterion is
    {                                                                // satisfied
      (*iteration_count_)++;                                         // increment iteration counter

      dolfin::symmetric_assemble(*matrix_, *matrixbc_, 
                                 *bilinear_, (*system_).dirichletbcs(), 
                                 NULL, NULL, NULL, 
                                 false);                             // assemble bilinear form
      dolfin::assemble(*rhs_, *linear_, false);                      // assemble linear form
      for(std::vector< const dolfin::DirichletBC* >::const_iterator  // loop over the collected vector of system bcs
                        bc = (*system_).dirichletbcs_begin(); 
                        bc != (*system_).dirichletbcs_end(); bc++)
      {
        (*(*bc)).apply(*rhs_);                                       // apply the bcs to the rhs 
      }
      (*matrixbc_).mult(*rhs_, *rhsbc_);
      *rhs_ -= *rhsbc_;

      for(std::vector<ReferencePoints_ptr>::const_iterator p =       // loop over the collected vector of system reference points
                                      (*system_).points_begin(); 
                                  p != (*system_).points_end(); p++)
      {
        (*(*p)).apply(*matrix_, *rhs_);                              // apply the reference points to the matrix and rhs
      }

      if(ident_zeros_)
      {
        (*matrix_).ident_zeros();
      }

      if (bilinearpc_)                                               // if there's a pc associated
      {
        assert(matrixpc_);
        dolfin::symmetric_assemble(*matrixpc_, *matrixbc_, 
                                   *bilinearpc_, (*system_).dirichletbcs(), 
                                   NULL, NULL, NULL, 
                                   false);                           // assemble the pc

        for(std::vector<ReferencePoints_ptr>::const_iterator p =     // loop over the collected vector of system reference points
                                        (*system_).points_begin(); 
                                    p != (*system_).points_end(); p++)
        {
          (*(*p)).apply(*matrixpc_);                                 // apply the reference points to the pc matrix
        }

        if(ident_zeros_pc_)
        {
          (*matrixpc_).ident_zeros();
        }

        perr = KSPSetOperators(ksp_, *(*matrix_).mat(),              // set the ksp operators with two matrices
                                     *(*matrixpc_).mat(), 
                                     SAME_NONZERO_PATTERN); 
        CHKERRV(perr);
      }
      else
      {
        perr = KSPSetOperators(ksp_, *(*matrix_).mat(),              // set the ksp operators with the same matrices
                                      *(*matrix_).mat(), 
                                        SAME_NONZERO_PATTERN); 
        CHKERRV(perr);
      }

      if (monitor_norms())
      {
        PetscReal norm;

        perr = VecNorm(*(*rhs_).vec(),NORM_2,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: 2-norm rhs = %f", norm);

        perr = VecNorm(*(*rhs_).vec(),NORM_INFINITY,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: inf-norm rhs = %f", norm);

        perr = VecNorm(*(*work_).vec(),NORM_2,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: 2-norm work = %f", norm);

        perr = VecNorm(*(*work_).vec(),NORM_INFINITY,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: inf-norm work = %f", norm);

        perr = MatNorm(*(*matrix_).mat(),NORM_FROBENIUS,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: Frobenius norm matrix = %f", norm);

        perr = MatNorm(*(*matrix_).mat(),NORM_INFINITY,&norm); CHKERRV(perr);
        dolfin::log(dolfin::get_log_level(), "Picard: inf-norm matrix = %f", norm);

        if (bilinearpc_)
        {
          perr = MatNorm(*(*matrixpc_).mat(),NORM_FROBENIUS,&norm); CHKERRV(perr);
          dolfin::log(dolfin::get_log_level(), "Picard: Frobenius norm matrix pc = %f", norm);

          perr = MatNorm(*(*matrixpc_).mat(),NORM_INFINITY,&norm); CHKERRV(perr);
          dolfin::log(dolfin::get_log_level(), "Picard: inf-norm matrix pc = %f", norm);
        }
      }

      perr = KSPSetUp(ksp_); CHKERRV(perr);                          // set up the ksp

      *work_ = (*(*(*system_).iteratedfunction()).vector());         // set the work vector to the iterated function
      perr = KSPSolve(ksp_, *(*rhs_).vec(), *(*work_).vec());        // perform a linear solve
      CHKERRV(perr);
      ksp_check_convergence_(ksp_);
      (*(*(*system_).iteratedfunction()).vector()) = *work_;         // update the iterated function with the work vector

      assert(residual_);
      dolfin::assemble(*res_, *residual_, false);                    // assemble the residual
      for(std::vector< const dolfin::DirichletBC* >::const_iterator bc = 
                             (*system_).dirichletbcs_begin(); 
                             bc != (*system_).dirichletbcs_end(); bc++)
      {                                                              // apply bcs to residual (should we do this?!)
        (*(*bc)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
      }
      for(std::vector<ReferencePoints_ptr>::const_iterator p = 
                                      (*system_).points_begin(); 
                                  p != (*system_).points_end(); p++)
      {                                                              // apply reference points to residual (should we do this?!)
        (*(*p)).apply(*res_, (*(*(*system_).iteratedfunction()).vector()));
      }

      aerror = (*res_).norm("l2");                                   // work out absolute error
      rerror = aerror/aerror0;                                       // and relative error
      dolfin::info("  %u Picard Residual Norm (absolute, relative) = %g, %g\n", 
                          iteration_count(), aerror, rerror);
                                                                     // and decide to loop or not...

      if(convfile_)
      {
        *(*(*system()).residualfunction()).vector() = (*boost::dynamic_pointer_cast< dolfin::GenericVector >(residual_vector()));
        if (convfile_)
        {
          (*convfile_).write_data();
        }
      }

    }

    if (iteration_count() == maxits_ && rerror > rtol_ && aerror > atol_)
    {
      dolfin::log(dolfin::WARNING, "it = %d, maxits_ = %d", iteration_count(), maxits_);
      dolfin::log(dolfin::WARNING, "rerror = %f, rtol_ = %f", rerror, rtol_);
      dolfin::log(dolfin::WARNING, "aerror = %f, atol_ = %f", aerror, atol_);
      if (ignore_failures_)
      {
        dolfin::log(dolfin::WARNING, "Picard iterations failed to converge, ignoring.");
      }
      else
      {
        dolfin::log(dolfin::ERROR, "Picard iterations failed to converge, sending sig int.");
        (*SignalHandler::instance()).dispatcher(SIGINT);
      }
    }

    (*(*(*system_).function()).vector()) =                              // update the function values with the iterated values
                      (*(*(*system_).iteratedfunction()).vector());

  }
  else                                                               // don't know what solver type this is
  {
    dolfin::error("Unknown solver type.");
  }

}

//*******************************************************************|************************************************************//
// assemble all linear forms (this includes initializing the vectors if necessary)
//*******************************************************************|************************************************************//
void SolverBucket::assemble_linearforms()
{
  assert(linear_);
  dolfin::assemble(*rhs_, *linear_, false);                          // and assemble it

  if(residual_)                                                      // do we have a residual_ form?
  {                                                                  // yes...
    dolfin::assemble(*res_, *residual_, false);                      // and assemble it
  }
}

//*******************************************************************|************************************************************//
// assemble all bilinear forms (this includes initializing the matrices if necessary)
//*******************************************************************|************************************************************//
void SolverBucket::assemble_bilinearforms()
{
  assert(bilinear_);
  dolfin::symmetric_assemble(*matrix_, *matrixbc_, 
                             *bilinear_, (*system_).dirichletbcs(), 
                             NULL, NULL, NULL, 
                             false, false, false);                   // assemble bilinear form
  Assembler::add_zeros_diagonal(*matrix_);                           // add zeros to the diagonal to ensure they remain in sparsity
  (*matrix_).apply("add");
  Assembler::add_zeros_diagonal(*matrixbc_);                         // add zeros to the diagonal to ensure they remain in sparsity
  (*matrixbc_).apply("add");

  if(bilinearpc_)                                                    // do we have a pc form?
  {
    dolfin::symmetric_assemble(*matrixpc_, *matrixbc_, 
                               *bilinear_, (*system_).dirichletbcs(), 
                               NULL, NULL, NULL, 
                               false, false, false);                   // assemble bilinear form
    Assembler::add_zeros_diagonal(*matrixpc_);                       // add zeros to the diagonal to ensure they remain in sparsity
    (*matrixpc_).apply("add");
    Assembler::add_zeros_diagonal(*matrixbc_);                         // add zeros to the diagonal to ensure they remain in sparsity
    (*matrixbc_).apply("add");

  }
}

//*******************************************************************|************************************************************//
// loop over the forms in this solver bucket and attach the coefficients they request using the parent bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::attach_form_coeffs()
{
  (*(*system_).bucket()).attach_coeffs(forms_begin(), forms_end());
}

//*******************************************************************|************************************************************//
// make a partial copy of the provided solver bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SolverBucket::copy_diagnostics(SolverBucket_ptr &solver, SystemBucket_ptr &system) const
{

  if(!solver)
  {
    solver.reset( new SolverBucket(&(*system)) );
  }

  (*solver).iteration_count_ = iteration_count_;
  (*solver).name_ = name_;
  (*solver).type_ = type_;   
  (*solver).copy_ = true;                                            // this is done to ensure that the petsc destroy routines
                                                                     // are not called by the destructor

}

//*******************************************************************|************************************************************//
// initialize any diagnostic output from the solver
//*******************************************************************|************************************************************//
void SolverBucket::initialize_diagnostics() const                    // doesn't allocate anything so can be const
{
  if (convfile_)
  {
    (*convfile_).write_header((*(*system()).bucket()));
  }
  if (kspconvfile_)
  {
    (*kspconvfile_).write_header((*(*system()).bucket()));
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
// return a pointer to the convergence file
//*******************************************************************|************************************************************//
const ConvergenceFile_ptr SolverBucket::convergence_file() const
{
  return convfile_;
}

//*******************************************************************|************************************************************//
// return a pointer to the ksp convergence file
//*******************************************************************|************************************************************//
const KSPConvergenceFile_ptr SolverBucket::ksp_convergence_file() const
{
  return kspconvfile_;
}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps
//*******************************************************************|************************************************************//
void SolverBucket::register_form(Form_ptr form, const std::string &name)
{
  Form_it f_it = forms_.find(name);                                  // check if this name already exists
  if (f_it != forms_.end())
  {
    dolfin::error("Form named \"%s\" already exists in solver.",     // if it does, issue an error
                                                  name.c_str());
  }
  else
  {
    forms_[name] = form;                                             // if not, register the form in the maps
  }
}

//*******************************************************************|************************************************************//
// return a boolean indicating if the solver bucket contains a form with the given name
//*******************************************************************|************************************************************//
bool SolverBucket::contains_form(const std::string &name)                   
{
  Form_it f_it = forms_.find(name);
  return f_it != forms_.end();
}

//*******************************************************************|************************************************************//
// return a (boost shared) pointer to a form from the solver bucket data maps
//*******************************************************************|************************************************************//
Form_ptr SolverBucket::fetch_form(const std::string &name)
{
  Form_it f_it = forms_.find(name);                                  // check if this name already exists
  if (f_it == forms_.end())
  {
    dolfin::error("Form named \"%s\" does not exist in solver.",     // if it doesn't, issue an error
                                                    name.c_str());
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
  return forms_.begin();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the beginning of the forms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::forms_begin() const
{
  return forms_.begin();
}

//*******************************************************************|************************************************************//
// return an iterator to the end of the forms_ map
//*******************************************************************|************************************************************//
Form_it SolverBucket::forms_end()
{
  return forms_.end();
}

//*******************************************************************|************************************************************//
// return a constant iterator to the end of the forms_ map
//*******************************************************************|************************************************************//
Form_const_it SolverBucket::forms_end() const
{
  return forms_.end();
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the solver bucket
//*******************************************************************|************************************************************//
const std::string SolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the forms in the solver bucket
//*******************************************************************|************************************************************//
const std::string SolverBucket::forms_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = forms_.begin(); f_it != forms_.end(); f_it++ )
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

  dolfin::log(dolfin::INFO, "Convergence for %s::%s", 
                          (*system_).name().c_str(), name().c_str());

  SNESConvergedReason snesreason;                                    // check what the convergence reason was
  PetscInt snesiterations;
  PetscInt sneslsiterations;
//  const char **snesprefix;
//  perr = SNESGetOptionsPrefix(snes_, snesprefix); CHKERRV(perr);   // FIXME: segfaults!
  perr = SNESGetConvergedReason(snes_, &snesreason); CHKERRV(perr);     
  perr = SNESGetIterationNumber(snes_, &snesiterations); CHKERRV(perr);
  perr = SNESGetLinearSolveIterations(snes_, &sneslsiterations);  CHKERRV(perr);
  dolfin::log(dolfin::INFO, "SNESConvergedReason %d", snesreason);
  dolfin::log(dolfin::INFO, "SNES n/o iterations %d", 
                              snesiterations);
  dolfin::log(dolfin::INFO, "SNES n/o linear solver iterations %d", 
                              sneslsiterations);
  if (snesreason<0)
  {
    if (ignore_failures_)
    {
      dolfin::log(dolfin::WARNING, "SNESConvergedReason <= 0, ignoring.");
    }
    else
    {
      dolfin::log(dolfin::ERROR, "SNESConvergedReason <= 0, sending sig int.");
      (*SignalHandler::instance()).dispatcher(SIGINT);
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
    dolfin::log(dolfin::INFO, "Convergence for %s::%s", 
                          (*system_).name().c_str(), name().c_str());
  }

  KSPConvergedReason kspreason;                                      // check what the convergence reason was
  PetscInt kspiterations;
  //const char **kspprefix;
  //perr = KSPGetOptionsPrefix(ksp, kspprefix); CHKERRV(perr);       // FIXME: segfaults!
  perr = KSPGetConvergedReason(ksp, &kspreason); CHKERRV(perr);     
  perr = KSPGetIterationNumber(ksp, &kspiterations); CHKERRV(perr);     
  dolfin::log(dolfin::INFO, "%sKSPConvergedReason %d", 
                              indentation.c_str(), kspreason);
  dolfin::log(dolfin::INFO, "%sKSP n/o iterations %d", 
                              indentation.c_str(), kspiterations);
  if (kspreason<0)
  {
    if (ignore_failures_)
    {
      dolfin::log(dolfin::WARNING, "KSPConvergedReason <= 0, ignoring.");
    }
    else
    {
      dolfin::log(dolfin::ERROR, "KSPConvergedReason <= 0, sending sig int.");
      (*SignalHandler::instance()).dispatcher(SIGINT);
    }
  }


  indent++;

  PC pc;
  const PCType pctype;
  perr = KSPGetPC(ksp, &pc); CHKERRV(perr);
  perr = PCGetType(pc, &pctype); CHKERRV(perr);

  if ((std::string)pctype=="ksp")
  {
    KSP subksp;                                                      // get the subksp from this pc
    perr = PCKSPGetKSP(pc, &subksp); CHKERRV(perr);
    ksp_check_convergence_(subksp, indent);
  }
  else if ((std::string)pctype=="fieldsplit")
  {
    KSP *subksps;                                                    // get the fieldsplit subksps
    PetscInt nsubksps;
    perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); 
    CHKERRV(perr); 
    for (PetscInt i = 0; i < nsubksps; i++)
    {
      ksp_check_convergence_(subksps[i], indent);
    }
  }
  
}

//*******************************************************************|************************************************************//
// empty the data structures in the solver bucket
//*******************************************************************|************************************************************//
void SolverBucket::empty_()
{
  forms_.clear();
}

