
#include "BoostTypes.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemSolversWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "petscsnes.h"

using namespace buckettools;

//*******************************************************************|************************************************************//
// specific constructor
//*******************************************************************|************************************************************//
SpudSolverBucket::SpudSolverBucket(const std::string &optionpath, 
                                            SystemBucket* system) : 
                                            optionpath_(optionpath), 
                                            SolverBucket(system)
{
                                                                     // do nothing
}

//*******************************************************************|************************************************************//
// default destructor
//*******************************************************************|************************************************************//
SpudSolverBucket::~SpudSolverBucket()
{
  empty_();                                                          // empty the data in the derived class
}

//*******************************************************************|************************************************************//
// fill the solver bucket data structures assuming the buckettools schema
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  fill_base_();                                                      // fill the base solver data: type, name etc.

  fill_forms_();                                                     // fill the forms data

  std::stringstream prefix;                                          // prefix buffer
  prefix.str(""); prefix << (*system_).name() << "_" << name() << "_";

  if (type()=="SNES")                                                // if this is a snes solver.  FIXME: switch to enum check
  {

    perr = SNESCreate(PETSC_COMM_WORLD, &snes_); CHKERRV(perr);      // create the petsc snes object

    perr = SNESSetOptionsPrefix(snes_, prefix.str().c_str());        // set its petsc options name prefix to SystemName_SolverName
    CHKERRV(perr);

    perr = SNESSetFromOptions(snes_); CHKERRV(perr);                 // set-up snes from options (we do this first to ensure that
                                                                     // any duplicated options from the options file overwrite the
                                                                     // command line)

    std::string snestype;
    buffer.str(""); buffer << optionpath() << "/type/snes_type/name";
    serr = Spud::get_option(buffer.str(), snestype);                 // set the snes type... ls is most common
    spud_err(buffer.str(), serr);
    perr = SNESSetType(snes_, snestype.c_str()); CHKERRV(perr); 

    buffer.str(""); buffer << optionpath() 
                                        << "/type/monitors/residual";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorDefault,               // set a snes residual monitor
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath() 
                                  << "/type/monitors/solution_graph";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorSolution,              // set a snes solution monitor (graph)
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    perr = SNESSetTolerances(snes_, atol_, rtol_, stol_, maxits_,    // from the data we collected in the base data fill set-up the
                                                            maxfes_);// snes tolerances
    CHKERRV(perr);

    perr = SNESGetKSP(snes_, &ksp_); CHKERRV(perr);                  // we always have at least one ksp so use the solverbucket ksp
                                                                     // to start setting up the ksp inside the snes
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // the ksp solver path
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // can then be used to fill the ksp data

    perr = SNESView(snes_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr); // turn on snesview so we get some debugging info

  }
  else if (type()=="Picard")                                         // if this is a picard solver
  {

    perr = KSPCreate(PETSC_COMM_WORLD, &ksp_); CHKERRV(perr);        // create a ksp object from the variable in the solverbucket

    buffer.str(""); buffer << optionpath() << "/type/linear_solver"; // figure out the linear solver optionspath
    fill_ksp_(buffer.str(), ksp_, prefix.str());                     // fill the ksp data

    perr = KSPView(ksp_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);   // turn on kspview so we get some debugging info

  }
  else                                                               // don't know how we got here
  {
    dolfin::error("Unknown solver type.");
  }

}

//*******************************************************************|************************************************************//
// make a partial copy of the provided solver bucket with the data necessary for writing the diagnostics file(s)
//*******************************************************************|************************************************************//
void SpudSolverBucket::copy_diagnostics(SolverBucket_ptr &solver, SystemBucket_ptr &system) const
{

  if(!solver)
  {
    solver.reset( new SpudSolverBucket(optionpath_, &(*system)) );
  }

  SolverBucket::copy_diagnostics(solver, system);

  (*boost::dynamic_pointer_cast< SpudSolverBucket >(solver)).form_optionpaths_ = form_optionpaths_;

}

//*******************************************************************|************************************************************//
// register a (boost shared) pointer to a form in the solver bucket data maps (with an optionpath as well)
//*******************************************************************|************************************************************//
void SpudSolverBucket::register_form(Form_ptr form, 
                                      const std::string &name, 
                                      const std::string &optionpath)
{
  Form_it f_it = forms_.find(name);                                  // check if name exists
  if (f_it != forms_.end())
  {
    dolfin::error("Form named \"%s\" already exists in function.",   // if it does, issue an error
                                                    name.c_str());
  }
  else
  {
    forms_[name]            = form;                                  // if not, insert form pointer into data map
    form_optionpaths_[name] = optionpath;                            // and do the same for its optionpath
  }
}

//*******************************************************************|************************************************************//
// return a string describing the contents of the solver bucket
//*******************************************************************|************************************************************//
const std::string SpudSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << " (" << 
                                    optionpath() << ")" << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

//*******************************************************************|************************************************************//
// return a string describing the forms in the solver bucket
//*******************************************************************|************************************************************//
const std::string SpudSolverBucket::forms_str(const int &indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = form_optionpaths_.begin(); 
                            s_it != form_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << 
                                (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

//*******************************************************************|************************************************************//
// fill the solver bucket base data assuming the buckettools schema (common for all solver types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_base_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code

  buffer.str(""); buffer << optionpath() << "/name";                 // solver name
  serr = Spud::get_option(buffer.str(), name_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/name";            // solver type (as string)
  serr = Spud::get_option(buffer.str(), type_);                      // FIXME: add conversion to enum here
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/relative_error";  // relative nonlinear error
  serr = Spud::get_option(buffer.str(), rtol_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/absolute_error";  // absolute nonlinear error
  serr = Spud::get_option(buffer.str(), atol_, 1.e-50); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/solution_error";  // nonlinear solution error (only applies to snes solver types)
  serr = Spud::get_option(buffer.str(), stol_, 1.e-8); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/max_iterations";  // maximum number of nonlinear iterations
  serr = Spud::get_option(buffer.str(), maxits_); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/min_iterations";  // minimum number of nonlinear iterations (only applies to
  serr = Spud::get_option(buffer.str(), minits_, 0);                 // picard solver types)
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() <<                          // maximum number of residual evaluations (only applies to snes
                                  "/type/max_function_evaluations";  // solver types)
  serr = Spud::get_option(buffer.str(), maxfes_, 10000); 
  spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << 
                                "/type/ignore_all_solver_failures";
  ignore_failures_ = Spud::have_option(buffer.str());

}

//*******************************************************************|************************************************************//
// fill the forms in solver bucket assuming the buckettools schema (common for all solver types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_forms_()
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
   
  buffer.str(""); buffer << optionpath() << "/type/form";           
  int nforms = Spud::option_count(buffer.str());                     // find out how many forms we have
  for (uint i = 0; i < nforms; i++)                                  // loop over them 
  {
    buffer.str(""); buffer << optionpath() << "/type/form[" << 
                                                          i << "]";
    std::string formoptionpath = buffer.str();

    std::string formname;
    buffer.str(""); buffer << formoptionpath << "/name";             // get the form name
    serr = Spud::get_option(buffer.str(), formname); 
    spud_err(buffer.str(), serr);

    Form_ptr form = ufc_fetch_form((*system_).name(), name(),        // get the form from the ufc using the system functionspace
                                        type(), 
                                        formname, 
                                        (*system_).functionspace());
    register_form(form, formname, formoptionpath);

                                                                     // at this stage we cannot attach any coefficients to this
                                                                     // form because we do not necessarily have them all
                                                                     // initialized yet so for the time being let's just grab any
                                                                     // functionspaces for the coefficients that we can find...
    uint ncoeff = (*form).num_coefficients();                        // how many coefficients does this form require?
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*form).coefficient_name(i);           // what is the (possibly derived) ufl symbol for this
                                                                     // coefficient
      if ((*(*system_).bucket()).contains_baseuflsymbol(uflsymbol))  // a base ufl symbol was only inserted into the parent bucket's
      {                                                              // if this is a coefficient function so we use this as an
                                                                     // indicator or whether we need to grab the functionspace or
                                                                     // not...
        std::string baseuflsymbol =                                  // what is the base ufl symbol?
              (*(*system_).bucket()).fetch_baseuflsymbol(uflsymbol); // have we already registered a functionspace for this base ufl
                                                                     // symbol?
        if (!(*(*system_).bucket()).contains_coefficientspace(baseuflsymbol))
        {                                                            // no...
          FunctionSpace_ptr coefficientspace;
          coefficientspace = ufc_fetch_coefficientspace_from_solver( // take a pointer to the functionspace from the ufc
                                        (*system_).name(), name(), 
                                        baseuflsymbol, 
                                        (*system_).mesh());
          (*(*system_).bucket()).register_coefficientspace(          // and register it in the parent bucket's map
                                        coefficientspace, 
                                        baseuflsymbol);
        }
      }

    }

  }

                                                                     // depending on the type of solver we assign certain forms to
                                                                     // hardcoded names in the solver bucket (basically it depends
                                                                     // which are linear or bilinear)
  ident_zeros_ = false;
  ident_zeros_pc_ = false;

  if (type()=="SNES")                                                // snes solver type...
  {
    linear_      = fetch_form("Residual");
    bilinear_    = fetch_form("Jacobian");

    buffer.str(""); buffer << optionpath() << "/type::SNES/form::Jacobian/ident_zeros";
    ident_zeros_ = Spud::have_option(buffer.str());

    if (contains_form("JacobianPC"))                                 // is there a pc form?
    {
      bilinearpc_ = fetch_form("JacobianPC");                        // yes but

      buffer.str(""); buffer << optionpath() << "/type::SNES/form::JacobianPC/ident_zeros";
      ident_zeros_pc_ = Spud::have_option(buffer.str());

    }                                                                // otherwise bilinearpc_ is null (indicates self pcing)
                                                                     // residual_ is always a null pointer for snes

  }
  else if (type()=="Picard")                                         // picard solver type...
  {
    linear_      = fetch_form("Linear");
    bilinear_    = fetch_form("Bilinear");

    buffer.str(""); buffer << optionpath() << "/type::Picard/form::Bilinear/ident_zeros";
    ident_zeros_ = Spud::have_option(buffer.str());

    if (contains_form("BilinearPC"))                                 // is there a pc form?
    {
      bilinearpc_ = fetch_form("BilinearPC");                        // yes but

      buffer.str(""); buffer << optionpath() << "/type::Picard/form::BilinearPC/ident_zeros";
      ident_zeros_pc_ = Spud::have_option(buffer.str());

    }                                                                // otherwise bilinearpc_ is null (indicates self pcing)
    residual_   = fetch_form("Residual");

  }
  else                                                               // unknown solver type
  {
    dolfin::error("Unknown solver type.");
  }
     
}

//*******************************************************************|************************************************************//
// fill a ksp object from the options tree (not necessarily the main solver bucket ksp_ object as this routine may be called
// recursively for ksp and fieldsplit pc types)
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                 const std::string prefix,
                                 const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  perr = KSPSetOptionsPrefix(ksp, prefix.c_str());                   // set the ksp options prefix
  CHKERRV(perr);

  std::string iterative_method;
  buffer.str(""); buffer << optionpath << "/iterative_method/name";  // iterative method (gmres, fgmres, cg etc.)
  serr = Spud::get_option(buffer.str(), iterative_method); 
  spud_err(buffer.str(), serr);

  perr = KSPSetType(ksp, iterative_method.c_str()); CHKERRV(perr);   // set the ksp type to the iterative method

  if(iterative_method != "preonly")                                  // tolerances (and monitors) only apply to iterative methods
  {
    PetscReal rtol;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/relative_error";
    serr = Spud::get_option(buffer.str(), rtol);                     // relative error
    spud_err(buffer.str(), serr);

    PetscReal atol;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/absolute_error";
    serr = Spud::get_option(buffer.str(), atol, 1.e-50);             // absolute error
    spud_err(buffer.str(), serr);

    PetscReal dtol;
    buffer.str(""); buffer << optionpath << 
                                "/iterative_method/divergence_error";
    serr = Spud::get_option(buffer.str(), dtol, 10000.0);            // divergence error (tolerance?) - preconditioned solution must
    spud_err(buffer.str(), serr);                                    // diverge this much to be considered divergent

    PetscInt maxits;
    buffer.str(""); buffer << optionpath << 
                                  "/iterative_method/max_iterations";// maximum number of linear solver iterations
    serr = Spud::get_option(buffer.str(), maxits); 
    spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << 
                "/iterative_method/monitors/preconditioned_residual";// monitor the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorDefault, 
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << 
                          "/iterative_method/monitors/true_residual";// monitor the true residual (more expensive)
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, 
                                            PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << 
          "/iterative_method/monitors/preconditioned_residual_graph";// plot a graph of the preconditioned residual
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorLG, 
                                             PETSC_NULL, PETSC_NULL); 
      CHKERRV(perr);
    }

    perr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  }

  buffer.str(""); buffer << optionpath << "/remove_null_space";      // removing a (or multiple) null space(s)
  if (Spud::have_option(buffer.str()))
  {
    buffer.str(""); buffer << optionpath << 
                                "/remove_null_space/null_space";
    int nnulls = Spud::option_count(buffer.str());                   // how many null spaces?

    uint kspsize = 0;
    if(parent_indices)
    {                                                                // if parent_indices is associated then this is a null space
       kspsize = (*parent_indices).size();                           // of a subksp so the kspsize is not the whole thing
    }                                                                // FIXME: broken in parallel!
    else
    {                                                                 
      kspsize = (*(*(*system_).function()).vector()).size();         // otherwise, this is quite easy - just the size of the parent
    }                                                                // system function

    std::vector< PETScVector_ptr > nullvecs;                         // collect the null space vectors here (so we maintain a reference)
    Vec vecs[nnulls];                                                // and here (for the petsc interface)
    std::vector<uint> prev_indices;

    for (uint i = 0; i<nnulls; i++)                                  // loop over the nullspaces
    {
      IS is;
      std::vector<uint> indices;                                     // record the indices from each iteration so that overlapping
                                                                     // null spaces may be avoided

      buffer.str(""); buffer << optionpath <<                        // optionpath of the nullspace
                        "/remove_null_space/null_space[" << i << "]";
      if (i==0)
      {
        fill_is_by_field_(buffer.str(), is,                          // create an is based on this optionpath (consistent schema
                          indices, parent_indices, NULL);            // with fieldsplit description)
      }
      else
      {
        fill_is_by_field_(buffer.str(), is,                          // create an is based on this optionpath (consistent schema
                          indices, parent_indices, &prev_indices);   // with fieldsplit description)
      }
      prev_indices.clear();
      prev_indices = indices;

      PETScVector_ptr nullvec( new dolfin::PETScVector(kspsize) );   // create a null vector for this null space
    
      PetscInt size, localsize;                                      // get the local and global sizes of this IS
      ISGetSize(is, &size);
      ISGetLocalSize(is, &localsize);
      dolfin::PETScVector unitvec(localsize, "local");               // create a local vector of local size length 
      double unit = 1./(std::sqrt( (double) size));                  // figure out what value a unit vector would have
      unitvec = unit;                                                // and set all entries of the vector to it

      (*nullvec).zero();                                             // zero the null vector

      VecScatter scatter;                                            // create a petsc scatter object from an object with the same 
      perr = VecScatterCreate(*unitvec.vec(), PETSC_NULL,            // structure as the unit vector to one with the same structure
                              *(*nullvec).vec(), is, &scatter);      // as the null vector using the IS
      CHKERRV(perr);
      perr = VecScatterBegin(scatter, *unitvec.vec(),                // scatter from the unit vector to the null vector
                             *(*nullvec).vec(), INSERT_VALUES, 
                             SCATTER_FORWARD); 
      CHKERRV(perr);
      perr = VecScatterEnd(scatter, *unitvec.vec(), 
                           *(*nullvec).vec(), INSERT_VALUES, 
                           SCATTER_FORWARD); 
      CHKERRV(perr);
      perr = VecScatterDestroy(scatter); CHKERRV(perr);              // necessary or taken care of when object leaves scope?
      perr = ISDestroy(is); CHKERRV(perr);                           // necessary or taken care of when object leaves scope?
      (*nullvec).apply("");                                          // finish assembly of the null vector, just in case

      nullvecs.push_back(nullvec);                                   // keep the null vector in scope by grabbing a reference to it
      vecs[i] = *(*nullvec).vec();                                   // also collect it in a petsc compatible format (shouldn't take
                                                                     // reference though... hence line above, necessary?)

    }

    MatNullSpace SP;                                                 // create a set of nullspaces in a null space object
    perr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, nnulls, 
                                                        vecs, &SP); 
    CHKERRV(perr);
    perr = KSPSetNullSpace(ksp, SP); CHKERRV(perr);                  // attach it to the ksp
    perr = MatNullSpaceDestroy(SP); CHKERRV(perr);                   // destroy the null space object, necessary?
  }

  std::string preconditioner;
  buffer.str(""); buffer << optionpath << "/preconditioner/name";    // preconditioner type
  serr = Spud::get_option(buffer.str(), preconditioner); 
  spud_err(buffer.str(), serr);

  PC pc;
  perr = KSPGetPC(ksp, &pc); CHKERRV(perr);                          // get the pc from the ksp
  perr = PCSetType(pc, preconditioner.c_str()); CHKERRV(perr);       // set its type (read from options earlier)

  if (preconditioner=="ksp")                                         // if the pc is itself a ksp
  {
    buffer.str(""); buffer << optionpath << 
                                    "/preconditioner/linear_solver";
    KSP subksp;                                                      // create a subksp from this pc
    perr = PCKSPGetKSP(pc, &subksp); CHKERRV(perr);
    fill_ksp_(buffer.str(), subksp, prefix+"subksp_", parent_indices);// recursively fill the ksp data (i.e. go back to this routine)
  }
  else if (preconditioner=="fieldsplit")                             // if the pc is a fieldsplit
  {
    buffer.str(""); buffer << optionpath << "/preconditioner";
    fill_pc_fieldsplit_(buffer.str(), pc, prefix, parent_indices);   // fill the fieldsplit data (will end up back here again
                                                                     // eventually)
  }
  else if (preconditioner=="lu")                                     // if the pc is direct
  {
    std::string factorization_package;                               // we get to choose a factorization package
    buffer.str(""); buffer << optionpath << 
                        "/preconditioner/factorization_package/name";
    serr = Spud::get_option(buffer.str(), factorization_package); 
    spud_err(buffer.str(), serr);

    perr = PCFactorSetMatSolverPackage(pc, 
                                      factorization_package.c_str()); 
    CHKERRV(perr);

  }

}

//*******************************************************************|************************************************************//
// fill a pc object from the options tree assuming its a fieldsplit pc
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_pc_fieldsplit_(const std::string &optionpath, 
                                           PC &pc, const std::string prefix,
                                           const std::vector<uint>* parent_indices)
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  std::string ctype;
  buffer.str(""); buffer << optionpath << "/composite_type/name";    // composite type of fieldsplit (additive, multiplicative etc.)
  serr = Spud::get_option(buffer.str(), ctype);                      // sadly no string based interface provided to this (I think)
  spud_err(buffer.str(), serr);                                      // so hard code an if block

  if (ctype == "additive")                                           // additive fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "multiplicative")                                // multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "symmetric_multiplicative")                      // symmetric multiplicative fieldsplit
  {
    perr = PCFieldSplitSetType(pc, 
                          PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE); 
    CHKERRV(perr);
  }
  else if (ctype == "special")                                       // special fieldsplit (whatever that means!)
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SPECIAL); 
    CHKERRV(perr);
  }
  else if (ctype == "schur")                                         // schur fieldsplit
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR); 
    CHKERRV(perr);
  }
  else                                                               // unknown (to buckettools) fieldsplit composite type
  {
    dolfin::error("Unknown PCCompositeType.");
  }

  std::vector< std::vector<uint> > child_indices;                    // a vector of vectors to collect the child indices (the
                                                                     // subsets of the parent_indices vector (if associated) that
                                                                     // will themselves become the parent_indices on the next
                                                                     // recursion)

  buffer.str(""); buffer << optionpath << "/fieldsplit";
  int nsplits = Spud::option_count(buffer.str());                    // how many fieldsplits exist for this pc
  for (uint i = 0; i < nsplits; i++)                                 // loop over them all
  {
    buffer.str(""); buffer << optionpath << 
                                        "/fieldsplit[" << i << "]";
    std::vector<uint> indices;
    IS is;
    if (i==0)
    {
      fill_is_by_field_(buffer.str(), is,                            // setup an IS for each fieldsplit
                        indices, parent_indices,
                        NULL);
    }
    else
    {
      fill_is_by_field_(buffer.str(), is,                            // setup an IS for each fieldsplit
                        indices, parent_indices,
                        &(child_indices[i-1]));
    }
    perr = PCFieldSplitSetIS(pc, is); CHKERRV(perr);                 // set the fs using that IS
    perr = ISDestroy(is); CHKERRV(perr);                             // destroy the IS, necessary?
    child_indices.push_back(indices);                                // record the indices of the global vector that made it into
  }                                                                  // IS (and hence the fieldsplit)

  KSP *subksps;                                                      // setup the fieldsplit subksps
  PetscInt nsubksps;
  perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); 
  CHKERRV(perr); 

  assert(nsubksps==nsplits);

  for (uint i = 0; i < nsplits; i++)                                 // loop over the splits again
  {
    std::string fsname;
    buffer.str(""); buffer << optionpath << "/fieldsplit[" 
                                        << i << "]/name";
    serr = Spud::get_option(buffer.str(), fsname);
    spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << "/fieldsplit[" 
                                        << i << "]/linear_solver";
    fill_ksp_(buffer.str(), subksps[i], prefix+fsname+"_", 
                                                &child_indices[i]);  // recurse and fill in the data on each subksp but passing down
                                                                     // the child_indices as the new parent_indices
  }

}

//*******************************************************************|************************************************************//
// fill a petsc is object from the options tree (for fieldsplits and null spaces - must share common schema tree)
// IS's may be set up by field name, components of the field and regions of the domain (FIXME: also need surface ids)
// this leads to four combinations... fields, fields with components, fields with regions, fields with both components and regions
// additionally the resulting IS is checked for consistency with any parents in the tree and will fail if a full description is
// not given
//*******************************************************************|************************************************************//
void SpudSolverBucket::fill_is_by_field_(const std::string &optionpath, IS &is, 
                                         std::vector<uint> &child_indices, 
                                         const std::vector<uint>* parent_indices,
                                         const std::vector<uint>* sibling_indices)
{

  std::stringstream buffer;                                          // optionpath buffer
  Spud::OptionError serr;                                            // spud error code
  PetscErrorCode perr;                                               // petsc error code

  std::pair<uint, uint> ownership_range =                            // the parallel ownership range of the system functionspace
          (*(*(*system_).functionspace()).dofmap()).ownership_range();

  buffer.str(""); buffer << optionpath << "/field";                  // loop over the fields used to describe this IS
  int nfields = Spud::option_count(buffer.str());

  if (nfields==0)                                                    // if no fields have been specified...
  {
    boost::unordered_set<uint>  dof_set =                            // take the set of dof for this system
            (*(*(*system_).functionspace()).dofmap()).dofs();
    for (boost::unordered_set<uint>::const_iterator                  // loop over the dof in the set
                                    dof_it = dof_set.begin(); 
                                    dof_it != dof_set.end(); 
                                    dof_it++)
    {                                                                // and insert them into the child_indices vector
      if ((*dof_it >= ownership_range.first) &&                      // but first check that this process owns them
                            (*dof_it < ownership_range.second))      // (in parallel)
      {
        child_indices.push_back(*dof_it);
      }
    }
    
  }

  bool mixedsystem = (((*system_).fields_size())>1);

  for (uint i = 0; i < nfields; i++)                                 // loop over the fields that have been specified
  {                                                                  // (won't do anything in the case of nfields==0, where we've
    buffer.str(""); buffer << optionpath << "/field[" << i << "]";   //  already included all dofs)

    std::string fieldname;
    buffer.str(""); buffer << optionpath << "/field[" << i <<        // get the field name
                                                          "]/name";
    serr = Spud::get_option(buffer.str(), fieldname); 
    spud_err(buffer.str(), serr);
    
    int fieldindex = (*(*system_).fetch_field(fieldname)).index();   // using the name, get the field index

    buffer.str(""); buffer << optionpath << "/field[" << i <<        // are region id restrictions specified under this field?
                                                   "]/region_ids";
    if (Spud::have_option(buffer.str()))
    {                                                                // yes... 
      std::vector< int > region_ids;
      serr = Spud::get_option(buffer.str(), region_ids);             // get the region ids
      spud_err(buffer.str(), serr);

      Mesh_ptr mesh = (*system_).mesh();                             // get the mesh
      MeshFunction_uint_ptr cellidmeshfunction =                     // and the region id mesh function
                      (*mesh).domains().cell_domains(*mesh);

      boost::unordered_set<uint> dof_set;                            // set up an unordered set of dof

      buffer.str(""); buffer << optionpath << "/field[" << i <<      // are specific components of this field requested?
                                                    "]/components";
      if (Spud::have_option(buffer.str()))
      {                                                              // yes... **field+region+component**

        std::vector< int > components;
        serr = Spud::get_option(buffer.str(), components);           // get the requested components
        spud_err(buffer.str(), serr);

        std::vector<int>::iterator max_comp_it = 
              std::max_element(components.begin(), components.end());// check this component exists
        if (mixedsystem)
        {
          assert(*max_comp_it < (*(*(*(*system_).functionspace())[fieldindex]).element()).num_sub_elements());
        }
        else
        {
          assert(*max_comp_it < (*(*(*system_).functionspace()).element()).num_sub_elements());
        }
        
        for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)  // loop over the cells in the mesh
        {

          int cellid = (*cellidmeshfunction)[(*cell).index()];       // get the cell region id from the mesh function

          for (std::vector<int>::const_iterator id =                 // loop over the regions that have been requested
                                  region_ids.begin(); 
                                  id != region_ids.end(); id++)
          {
            if(cellid==*id)                                          // and check if this cell should be included
            {                                                        // yes...
              for (std::vector<int>::const_iterator comp =           // loop over the components that have been requested
                                      components.begin(); 
                                      comp != components.end(); 
                                      comp++)
              {
                std::vector<uint> dof_vec;
                if (mixedsystem)
                {
                  dof_vec =                                          // get the cell dof for the current component
                      (*(*(*(*(*system_).functionspace())[fieldindex])[*comp]).dofmap()).cell_dofs((*cell).index());
                }
                else
                {
                  dof_vec =                                          // get the cell dof for the current component
                      (*(*(*(*system_).functionspace())[*comp]).dofmap()).cell_dofs((*cell).index());
                }
                for (std::vector<uint>::const_iterator dof_it =      // loop over the cell dof for this component
                                              dof_vec.begin(); 
                                              dof_it < dof_vec.end(); 
                                              dof_it++)
                {
                  dof_set.insert(*dof_it);                           // insert each dof into the unordered set
                }                                                    // (i.e. if it's not there already)
              }
            }
          }
        }
      }
      else
      {                                                              // no, no components specified... **field+region**

        for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)  // loop over the cells in the mesh
        {
          int cellid = (*cellidmeshfunction)[(*cell).index()];       // get the cell region id from the mesh function

          for (std::vector<int>::const_iterator id =                 // loop over the region ids that have been requested
                                      region_ids.begin(); 
                                      id != region_ids.end(); id++)
          {
            if(cellid==*id)                                          // check if this cell should be included
            {                                                        // yes...
              std::vector<uint> dof_vec;
              if (mixedsystem)
              {
                dof_vec =                                            // get the cell dof (potentially for all components)
                  (*(*(*(*system_).functionspace())[fieldindex]).dofmap()).cell_dofs((*cell).index());
              }
              else
              {
                dof_vec =                                            // get the cell dof (potentially for all components)
                  (*(*(*system_).functionspace()).dofmap()).cell_dofs((*cell).index());
              }
              for (std::vector<uint>::const_iterator dof_it =        // loop over the cell dof
                                      dof_vec.begin(); 
                                      dof_it < dof_vec.end(); 
                                      dof_it++)
              {
                dof_set.insert(*dof_it);                             // and insert each one into the unordered set
              }                                                      // (i.e. if it hasn't been added already)
            }
          }
        }

      }
                                                                     // we now have a set of dofs for this field over all regions
                                                                     // and components requested...
                                                                     // put this into the vector of dofs for all fields (field
                                                                     // indices, unlike cell indices, should be unique so no need
                                                                     // for a set anymore)
      for (boost::unordered_set<uint>::const_iterator dof_it =      
                                            dof_set.begin(); 
                                            dof_it != dof_set.end(); // loop over all entries in the set 
                                            dof_it++)
      {                                                              // and insert into the child_indices vector
        if ((*dof_it >= ownership_range.first)                       // but only if the dof is within the ownership range in
                && (*dof_it < ownership_range.second))               // parallel
        {
          child_indices.push_back(*dof_it);
        }
      }
    }
    else
    {                                                                // no region ids specified, but there might still be components
      buffer.str(""); buffer << optionpath << "/field[" << i << 
                                                      "]/components";
      if (Spud::have_option(buffer.str()))
      {                                                              // yes, there are components specified... **field+components**
        std::vector< int > components;
        serr = Spud::get_option(buffer.str(), components);           // get the components
        spud_err(buffer.str(), serr);
        
        std::vector<int>::iterator max_comp_it =                  
             std::max_element(components.begin(), components.end()); // check this component exists
        if (mixedsystem)
        {
          assert(*max_comp_it < (*(*(*(*system_).functionspace())[fieldindex]).element()).num_sub_elements());
        }
        else
        {
          assert(*max_comp_it < (*(*(*system_).functionspace()).element()).num_sub_elements());
        }

        for (std::vector<int>::const_iterator comp =                 // loop over the requested components
                                      components.begin(); 
                                      comp != components.end(); 
                                      comp++)
        {
          boost::unordered_set<uint> dof_set;
          if (mixedsystem)
          {
            dof_set =                                                // take the set of dof for this component
                (*(*(*(*(*system_).functionspace())[fieldindex])[*comp]).dofmap()).dofs();
          }
          else
          {
            dof_set =                                                // take the set of dof for this component
                (*(*(*(*system_).functionspace())[*comp]).dofmap()).dofs();
          }
          for (boost::unordered_set<uint>::const_iterator            // loop over the dof in the set
                                        dof_it = dof_set.begin(); 
                                        dof_it != dof_set.end(); 
                                        dof_it++)
          {                                                          // and insert them into the child_indices vector
            if ((*dof_it >= ownership_range.first) &&                // but first check that this process owns them
                                (*dof_it < ownership_range.second))  // (in parallel)
            {
              child_indices.push_back(*dof_it);
            }
          }
        }
      }
      else
      {                                                              // no regions or components specified... **field**
        boost::unordered_set<uint> dof_set;
        if (mixedsystem)
        {
          dof_set =                                                  // take the set of dof for this field
                (*(*(*(*system_).functionspace())[fieldindex]).dofmap()).dofs();
        }
        else
        {
          dof_set =                                                  // take the set of dof for this field
                (*(*(*system_).functionspace()).dofmap()).dofs();
        }
        for (boost::unordered_set<uint>::const_iterator              // loop over the dof in the set
                                        dof_it = dof_set.begin(); 
                                        dof_it != dof_set.end(); 
                                        dof_it++)
        {                                                            // and insert them into the child_indices vector
          if ((*dof_it >= ownership_range.first) &&                  // but first check that this process owns them
                                (*dof_it < ownership_range.second))  // (in parallel)
          {
            child_indices.push_back(*dof_it);
          }
        }
      }
    }
  }
                                                                     // phew, that's over!
                                                                     // FIXME: tidy up the above mess!

  std::sort(child_indices.begin(), child_indices.end());             // sort the vector of child_indices

  if(sibling_indices)                                                // we have been passed a list of sibling indices...
  {                                                                  // we wish to remove from the child_indices any indices that
                                                                     // also occur in the sibling indices
    uint c_size = child_indices.size();
    uint c_ind = 0;
    std::vector<uint> tmp_child_indices;
    bool overlap = false;
    for(std::vector<uint>::const_iterator                            // loop over the sibling indices
                                   s_it = (*sibling_indices).begin();
                                   s_it != (*sibling_indices).end();
                                   s_it++)
    {
      while(child_indices[c_ind] != *s_it)                           // child_indices are sorted, so sibling_indices should be too
      {                                                              // search child_indices until the current sibling index is found
        tmp_child_indices.push_back(child_indices[c_ind]);           // include indices that aren't in the sibling
        c_ind++;
        if (c_ind == c_size)                                         // or we reach the end of the child_indices...
        {
          break;
        }
      }
      if (c_ind == c_size)                                           // we've reached the end of the child indices so nothing more
      {                                                              // to do
        break;
      }
      else                                                           // we haven't reached the end of the child indices but found
      {                                                              // a sibling index to ignore... give a warning
        overlap = true;
      }
      c_ind++;                                                       // indices shouldn't be repeated so incredment the child too
    }

    if(overlap)
    {                                                                // sibling indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices overlap with sibling fieldsplit, ignoring overlapping indices.");
    }
    child_indices.clear();
    child_indices = tmp_child_indices;
                                 
  }

  if(parent_indices)                                                 // we have been passed a list of parent indices... 
  {                                                                  // we wish to remove from the child_indices any indices that do
                                                                     // not occur in the parent indices 
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    uint p_reset = 0;
    std::vector<uint> tmp_child_indices;
    bool extra = false;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        c_it = child_indices.begin(); 
                                        c_it != child_indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // child_indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and prepare to throw a warning
          extra = true;
          break;
        }
      }
      if (p_ind == p_size)
      {
        p_ind = p_reset;
      }
      else
      {
        tmp_child_indices.push_back(*c_it);                          // include indices that are in the parent
        p_ind++;                                                     // indices shouldn't be repeated so increment the parent too
        p_reset = p_ind;                                             // this is where the next failed search should continue from
        if (p_ind == p_size)                                         // we've reached the end
        { 
          break;                                            
        }
      }
    } 

    if(extra)
    {                                                                // child indices were ignored... give a warning
      dolfin::log(dolfin::WARNING, 
                  "WARNING: IS indices not a subset of parent fieldsplit, ignoring extra indices.");
    }
    child_indices.clear();
    child_indices = tmp_child_indices;
                                 
  }

  PetscInt n=child_indices.size();                                   // setup a simpler structure for petsc
  assert(n>0);
  PetscInt *indices;
  PetscMalloc(n*sizeof(PetscInt), &indices);
 
  uint ind = 0;
  if(parent_indices)
  {                                                                  // we have been passed a list of parent indices... 
                                                                     // our child indices must be a  subset of this list and indexed
                                                                     // into it so let's do that now while we convert structures...
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    for (std::vector<uint>::const_iterator                           // loop over the child indices
                                        c_it = child_indices.begin(); 
                                        c_it != child_indices.end(); 
                                        c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)                      // child_indices is sorted, so parent_indices should be too...
      {                                                              // search parent_indices until the current child index is found
        p_ind++;
        if (p_ind == p_size)                                         // or we reach the end of the parent_indices...
        {                                                            // and throw an error
          dolfin::error("IS indices are not a subset of a parent fieldsplit, shouldn't happen here.");
        }
      }
      indices[ind] = p_ind;                                          // found the child index in the parent_indices so copy it into
                                                                     // the PetscInt array
      ind++;                                                         // increment the array index
      p_ind++;                                                       // indices shouldn't be repeated so increment the parent too
    } 
    assert(ind==n);                                                  // these should be equal
  }
  else
  {
    for (std::vector<uint>::const_iterator                           // loop over the child_indices
                                      ind_it = child_indices.begin(); 
                                      ind_it != child_indices.end(); 
                                      ind_it++)
    {
      indices[ind] = *ind_it;                                        // insert them into the PetscInt array
      ind++;                                                         // increment the array index
    }
    // these should be equal
    assert(ind==n);                                                  // these should be equal
  }

  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, &is);         // create the general index set based on the indices
  CHKERRV(perr);
  if (Spud::have_option(optionpath+"/monitors/view_index_set"))
  {
    buffer.str(""); buffer << optionpath << "/name";                 // IS Name
    std::string isname;
    serr = Spud::get_option(buffer.str(), isname);
    spud_err(buffer.str(), serr);
    
    dolfin::log(dolfin::INFO, "ISView: %s (%s)", 
                                isname.c_str(), optionpath.c_str());
    perr = ISView(is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);      // isview?
  }

  PetscFree(indices);                                                // free the PetscInt array of indices
   
}

//*******************************************************************|************************************************************//
// empty the data structures in the spudsolver bucket
//*******************************************************************|************************************************************//
void SpudSolverBucket::empty_()
{
  form_optionpaths_.clear();
  SolverBucket::empty_();
}

