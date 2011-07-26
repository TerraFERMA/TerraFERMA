
#include "BoostTypes.h"
#include "SolverBucket.h"
#include "SystemBucket.h"
#include <dolfin.h>
#include <string>

using namespace buckettools;

// Default constructor
SolverBucket::SolverBucket()
{
  // Do nothing
}

// Specific constructor
SolverBucket::SolverBucket(SystemBucket* system) : system_(system)
{
  // Do nothing
}

// Default destructor
SolverBucket::~SolverBucket()
{
  empty_();
}

// Register a form in the solver
void SolverBucket::register_form(Form_ptr form, std::string name)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it != forms_.end())
  {
    // if it does, issue an error
    dolfin::error("Form named \"%s\" already exists in solver.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    forms_[name] = form;
  }
}

// Do we have a form with a particular name?
bool SolverBucket::contains_form(std::string name)
{
  Form_it f_it = forms_.find(name);
  return f_it != forms_.end();
}

// Return a pointer to a form with the given name
Form_ptr SolverBucket::fetch_form(std::string name)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it == forms_.end())
  {
    // if it doesn't, issue an error
    dolfin::error("Form named \"%s\" does not exist in solver.", name.c_str());
  }
  else
  {
    // if not then return it
    return (*f_it).second;
  }
}

Form_it SolverBucket::forms_begin()
{
  return forms_.begin();
}

Form_const_it SolverBucket::forms_begin() const
{
  return forms_.begin();
}

Form_it SolverBucket::forms_end()
{
  return forms_.end();
}

Form_const_it SolverBucket::forms_end() const
{
  return forms_.end();
}

void SolverBucket::attach_form_coeffs()
{
  for (Form_it f_it = forms_begin(); f_it != forms_end(); f_it++)
  {
    // Loop over the functions requested by the form and hook up pointers
    uint ncoeff = (*(*f_it).second).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*(*f_it).second).coefficient_name(i);
      GenericFunction_ptr function = (*system_).fetch_uflsymbol(uflsymbol);

      (*(*f_it).second).set_coefficient(uflsymbol, function);
    }
  }

}

void SolverBucket::initialize_matrices()
{
  PetscErrorCode perr;

  if (type()=="SNES") 
  {
    ctx_.linear       = linear_;
    ctx_.bilinear     = bilinear_;
    ctx_.bilinearpc   = bilinearpc_;
    ctx_.bcs          = (*system_).bcs();
    ctx_.function     = (*system_).function();
    ctx_.oldfunction  = (*system_).oldfunction();
  }

  assemble_bilinearforms(true);
  assemble_linearforms(true);
  
  uint syssize = (*(*system_).function()).vector().size();
  work_.reset( new dolfin::PETScVector(syssize) ); 
  if(type()=="SNES")
  {
    assert(!res_);
    res_.reset( new dolfin::PETScVector(syssize) );

    perr = SNESSetFunction(snes_, *(*res_).vec(), FormFunction, (void *) &ctx_); CHKERRV(perr);

    if (bilinearpc_) 
    {
      assert(matrixpc_);
      perr = SNESSetJacobian(snes_, *(*matrix_).mat(), *(*matrixpc_).mat(), FormJacobian, (void *) &ctx_); CHKERRV(perr);
    }
    else 
    {
      perr = SNESSetJacobian(snes_, *(*matrix_).mat(), *(*matrix_).mat(), FormJacobian, (void *) &ctx_); CHKERRV(perr);
    }

    perr = SNESView(snes_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);

  }

}

void SolverBucket::solve()
{
  PetscErrorCode perr;

  if (type()=="SNES")
  {
    *work_ = (*(*system_).function()).vector();
    perr = SNESSolve(snes_, PETSC_NULL, *(*work_).vec());
    (*(*system_).function()).vector() = *work_;
  }
  else if (type()=="Picard")
  {

    uint it = 0;
    double aerror = (*res_).norm("l2");
    double aerror0 = aerror;
    double rerror = aerror/aerror0;
    dolfin::info("%u Error (absolute, relative) = %g, %g\n", it, aerror, rerror);

    while (it < minits_ || (it < maxits_ && rerror > rtol_ && aerror > atol_))
    {
      it++;

      dolfin::assemble(*matrix_, *bilinear_);
      dolfin::assemble(*rhs_, *linear_);
      for(std::vector<BoundaryCondition_ptr>::const_iterator bc = (*system_).bcs_begin(); bc != (*system_).bcs_end(); bc++)
      {
        (*(*bc)).apply(*matrix_, *rhs_);
      }

      if (bilinearpc_)
      {
        assert(matrixpc_);
        dolfin::assemble(*matrixpc_, *bilinearpc_);
        for(std::vector<BoundaryCondition_ptr>::const_iterator bc = (*system_).bcs_begin(); bc != (*system_).bcs_end(); bc++)
        {
          (*(*bc)).apply(*matrixpc_, *rhs_);
        }

        perr = KSPSetOperators(ksp_, *(*matrix_).mat(), *(*matrixpc_).mat(), SAME_NONZERO_PATTERN); CHKERRV(perr);
      }
      else
      {
        perr = KSPSetOperators(ksp_, *(*matrix_).mat(), *(*matrix_).mat(), SAME_NONZERO_PATTERN); CHKERRV(perr);
      }

      perr = KSPSetUp(ksp_); CHKERRV(perr);

      *work_ = (*(*system_).function()).vector();
      perr = KSPSolve(ksp_, *(*rhs_).vec(), *(*work_).vec());  CHKERRV(perr);
      (*(*system_).function()).vector() = *work_;

      assert(residual_);
      dolfin::assemble(*res_, *residual_);
      for(std::vector<BoundaryCondition_ptr>::const_iterator bc = (*system_).bcs_begin(); bc != (*system_).bcs_end(); bc++)
      {
        (*(*bc)).apply(*res_, (*(*system_).function()).vector());
      }

      aerror = (*res_).norm("l2");
      rerror = aerror/aerror0;
      dolfin::info("%u Error (absolute, relative) = %g, %g\n", it, aerror, rerror);

    }

  }
  else
  {
    dolfin::error("Unknown solver type.");
  }

}

void SolverBucket::assemble_bilinearforms(const bool &reset_tensor)
{

  assert(bilinear_);
  if(!matrix_)
  {
    matrix_.reset(new dolfin::PETScMatrix);
  }
  dolfin::assemble(*matrix_, *bilinear_, reset_tensor);
  if(bilinearpc_)
  {
    if(!matrixpc_)
    {
      matrixpc_.reset(new dolfin::PETScMatrix);
    }
    dolfin::assemble(*matrixpc_, *bilinearpc_, reset_tensor);
//    for(std::vector<BoundaryCondition_ptr>::const_iterator bc = (*system_).bcs_begin(); bc != (*system_).bcs_end(); bc++)
//    {
//      (*bc).apply(*res_, (*(*system_).function()))
//    }
  }

}

void SolverBucket::assemble_linearforms(const bool &reset_tensor)
{

  assert(linear_);
  if(!rhs_)
  {
    rhs_.reset(new dolfin::PETScVector);
  }
  dolfin::assemble(*rhs_, *linear_, reset_tensor);
  if(residual_)
  {
    if(!res_)
    {
      res_.reset(new dolfin::PETScVector);
    }
    dolfin::assemble(*res_, *residual_, reset_tensor);
//    for(std::vector<BoundaryCondition_ptr>::const_iterator bc = (*system_).bcs_begin(); bc != (*system_).bcs_end(); bc++)
//    {
//      (*bc).apply(*res_, (*(*system_).function()))
//    }
  }

}

// Return a string describing the contents of the function
const std::string SolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Return a string describing the contents of forms_
const std::string SolverBucket::forms_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  for ( Form_const_it f_it = forms_.begin(); f_it != forms_.end(); f_it++ )
  {
    s << indentation << "Form " << (*f_it).first  << std::endl;
  }
  return s.str();
}

// Empty the function
void SolverBucket::empty_()
{
  forms_.clear();
}

