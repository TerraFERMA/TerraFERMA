
#include "BoostTypes.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemsWrapper.h"
#include "SpudSystemBucket.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "petscsnes.h"

using namespace buckettools;

// Specific constructor
SpudSolverBucket::SpudSolverBucket(std::string optionpath, SystemBucket* system) : optionpath_(optionpath), SolverBucket(system)
{
  // Do nothing
}

// Default destructor (declared as virtual so will call base class destructor)
SpudSolverBucket::~SpudSolverBucket()
{
  PetscErrorCode perr;

  if(type()=="SNES")
  {
    perr = SNESDestroy(snes_); CHKERRV(perr);
  }

  if(type()=="Picard")
  {
    perr = KSPDestroy(ksp_); CHKERRV(perr);
  }

}

// Fill the function using spud and assuming a buckettools schema structure
void SpudSolverBucket::fill()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  base_fill_();

  forms_fill_();

  if (type()=="SNES")
  {
    perr = SNESCreate(PETSC_COMM_WORLD, &snes_); CHKERRV(perr);

    buffer.str(""); buffer << (*system_).name() << "_" << name() << "_";
    perr = SNESSetOptionsPrefix(snes_, buffer.str().c_str()); CHKERRV(perr);

    perr = SNESSetFromOptions(snes_); CHKERRV(perr);

    std::string snestype;
    buffer.str(""); buffer << optionpath() << "/type/snes_type/name";
    serr = Spud::get_option(buffer.str(), snestype); spud_err(buffer.str(), serr);

    perr = SNESSetType(snes_, snestype.c_str()); CHKERRV(perr);

    buffer.str(""); buffer << optionpath() << "/type/monitors/residual";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorDefault, PETSC_NULL, PETSC_NULL); CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath() << "/type/monitors/solution_graph";
    if (Spud::have_option(buffer.str()))
    {
      perr = SNESMonitorSet(snes_, SNESMonitorSolution, PETSC_NULL, PETSC_NULL); CHKERRV(perr);
    }

    perr = SNESSetTolerances(snes_, atol_, rtol_, stol_, maxits_, maxfes_); CHKERRV(perr);

    // we always have at least one ksp
    perr = SNESGetKSP(snes_, &ksp_); CHKERRV(perr);
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver";
    ksp_fill_(buffer.str(), ksp_);

    perr = SNESView(snes_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);
  }
  else if (type()=="Picard")
  {
    perr = KSPCreate(PETSC_COMM_WORLD, &ksp_); CHKERRV(perr);

    buffer.str(""); buffer << (*system_).name() << "_" << name() << "_";
    perr = KSPSetOptionsPrefix(ksp_, buffer.str().c_str()); CHKERRV(perr);

    buffer.str(""); buffer << optionpath() << "/type/linear_solver";
    ksp_fill_(buffer.str(), ksp_);

    perr = KSPView(ksp_, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);
  }
  else
  {
    dolfin::error("Unknown solver type.");
  }

}

void SpudSolverBucket::ksp_fill_(const std::string &optionpath, KSP &ksp, 
                                 const std::vector<uint>* parent_indices)
{
  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  std::string iterative_method;
  buffer.str(""); buffer << optionpath << "/iterative_method/name";
  serr = Spud::get_option(buffer.str(), iterative_method); spud_err(buffer.str(), serr);

  perr = KSPSetType(ksp, iterative_method.c_str()); CHKERRV(perr);

  if(iterative_method != "preonly")
  {
    PetscReal rtol;
    buffer.str(""); buffer << optionpath << "/iterative_method/relative_error";
    serr = Spud::get_option(buffer.str(), rtol); spud_err(buffer.str(), serr);

    PetscReal atol;
    buffer.str(""); buffer << optionpath << "/iterative_method/absolute_error";
    serr = Spud::get_option(buffer.str(), atol, 1.e-50); spud_err(buffer.str(), serr);

    PetscReal dtol;
    buffer.str(""); buffer << optionpath << "/iterative_method/divergence_error";
    serr = Spud::get_option(buffer.str(), dtol, 10000.0); spud_err(buffer.str(), serr);

    PetscInt maxits;
    buffer.str(""); buffer << optionpath << "/iterative_method/max_iterations";
    serr = Spud::get_option(buffer.str(), maxits); spud_err(buffer.str(), serr);

    buffer.str(""); buffer << optionpath << "/iterative_method/monitors/preconditioned_residual";
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL); CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << "/iterative_method/monitors/true_residual";
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL); CHKERRV(perr);
    }

    buffer.str(""); buffer << optionpath << "/iterative_method/monitors/preconditioned_residual_graph";
    if (Spud::have_option(buffer.str()))
    {
      perr = KSPMonitorSet(ksp, KSPMonitorLG, PETSC_NULL, PETSC_NULL); CHKERRV(perr);
    }

    perr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  }

  std::string preconditioner;
  buffer.str(""); buffer << optionpath << "/preconditioner/name";
  serr = Spud::get_option(buffer.str(), preconditioner); spud_err(buffer.str(), serr);

  PC pc;
  perr = KSPGetPC(ksp, &pc); CHKERRV(perr);
  perr = PCSetType(pc, preconditioner.c_str()); CHKERRV(perr);

  if (preconditioner=="ksp")
  {
    buffer.str(""); buffer << optionpath << "/preconditioner/linear_solver";
    KSP subksp;
    perr = PCKSPGetKSP(pc, &subksp); CHKERRV(perr);
    // recurse!
    ksp_fill_(buffer.str(), subksp, parent_indices);
  }
  else if (preconditioner=="fieldsplit")
  {
    buffer.str(""); buffer << optionpath << "/preconditioner";
    pc_fieldsplit_fill_(buffer.str(), pc, parent_indices);
  }  

}

void SpudSolverBucket::pc_fieldsplit_fill_(const std::string &optionpath, PC &pc, 
                                           const std::vector<uint>* parent_indices)
{

  std::stringstream buffer;
  PetscErrorCode perr;
  Spud::OptionError serr;

  std::string ctype;
  buffer.str(""); buffer << optionpath << "/composite_type/name";
  serr = Spud::get_option(buffer.str(), ctype); spud_err(buffer.str(), serr);
  if (ctype == "additive")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE); CHKERRV(perr);
  }
  else if (ctype == "multiplicative")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE); CHKERRV(perr);
  }
  else if (ctype == "symmetric_multiplicative")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE); CHKERRV(perr);
  }
  else if (ctype == "special")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SPECIAL); CHKERRV(perr);
  }
  else if (ctype == "schur")
  {
    perr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR); CHKERRV(perr);
  }
  else
  {
    dolfin::error("Unknown PCCompositeType.");
  }

  // create a vector of vectors to contain the
  // child indices (the subsets of the parent_indices vector 
  // (if associated))
  std::vector< std::vector<uint> > child_indices;

  buffer.str(""); buffer << optionpath << "/fieldsplit";
  int nsplits = Spud::option_count(buffer.str());
  for (uint i = 0; i < nsplits; i++)
  {
    buffer.str(""); buffer << optionpath << "/fieldsplit[" << i << "]";
    std::vector<uint> indices;
    pc_fieldsplit_by_field_fill_(buffer.str(), pc, indices, parent_indices);
    child_indices.push_back(indices);
  }

  KSP *subksps;
  PetscInt nsubksps;
  perr = PCFieldSplitGetSubKSP(pc, &nsubksps, &subksps); CHKERRV(perr); 

  assert(nsubksps==nsplits);

  for (uint i = 0; i < nsplits; i++)
  {
    buffer.str(""); buffer << optionpath << "/fieldsplit[" << i << "]/linear_solver";
    // recurse (but now child_indices[i] becomes the parent_indices for the next level!
    ksp_fill_(buffer.str(), subksps[i], &child_indices[i]);
  }

}

void SpudSolverBucket::pc_fieldsplit_by_field_fill_(const std::string &optionpath, PC &pc, 
                                                    std::vector<uint> &child_indices, 
                                                    const std::vector<uint>* parent_indices)
{

  // a buffer for strings
  std::stringstream buffer;
  // petsc error code holder
  PetscErrorCode perr;
  // spud error code holder
  Spud::OptionError serr;

  // the ownership range of this functionspace
  std::pair<uint, uint> ownership_range = (*(*system_).functionspace()).dofmap().ownership_range();

  buffer.str(""); buffer << optionpath << "/field";
  int nfields = Spud::option_count(buffer.str());
  // loop over the fields used in this fieldsplit
  for (uint i = 0; i < nfields; i++)
  {
    buffer.str(""); buffer << optionpath << "/field[" << i << "]";

    // get the fieldname
    std::string fieldname;
    buffer.str(""); buffer << optionpath << "/field[" << i << "]/name";
    serr = Spud::get_option(buffer.str(), fieldname); spud_err(buffer.str(), serr);
    
    // and from that the field index
    int fieldindex = (*(*system_).fetch_field(fieldname)).index();

    // do we specify region id limitations for this fieldsplit
    buffer.str(""); buffer << optionpath << "/field[" << i << "]/region_ids";
    if (Spud::have_option(buffer.str()))
    {
      // yes, get the region_ids
      std::vector< int > region_ids;
      serr = Spud::get_option(buffer.str(), region_ids); spud_err(buffer.str(), serr);

      // fetch the mesh and relevent mesh data
      Mesh_ptr mesh = (*system_).mesh();
      MeshFunction_uint_ptr cellidmeshfunction = (*mesh).data().mesh_function("CellIDs");

      // set up a set of dof
      boost::unordered_set<uint> dof_set;

      // do we specify components of this field?
      buffer.str(""); buffer << optionpath << "/field[" << i << "]/components";
      if (Spud::have_option(buffer.str()))
      {
        // yes, get the component list
        std::vector< int > components;
        serr = Spud::get_option(buffer.str(), components); spud_err(buffer.str(), serr);

        // check that the user hasn't requested a component that doesn't exist
        std::vector<int>::iterator max_comp_it = std::max_element(components.begin(), components.end());
        assert(*max_comp_it < (*(*(*system_).functionspace())[fieldindex]).element().num_sub_elements());
        
        // loop over the cells in the mesh
        for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)
        {
          // get the cell id from the mesh function
          int cellid = (*cellidmeshfunction)[(*cell).index()];
          // for each region_id we specified,
          for (std::vector<int>::const_iterator id = region_ids.begin(); id != region_ids.end(); id++)
          {
            // check if this cell should be included
            if(cellid==*id)
            {
              // if yes, then loop over the components we want included
              for (std::vector<int>::const_iterator comp = components.begin(); comp != components.end(); comp++)
              {
                // and get their dofs from the dofmap
                std::vector<uint> dof_vec = (*(*(*(*system_).functionspace())[fieldindex])[*comp]).dofmap().cell_dofs((*cell).index());
                for (std::vector<uint>::const_iterator dof_it = dof_vec.begin(); dof_it < dof_vec.end(); dof_it++)
                {
                  // insert them into the set of dof for this field
                  dof_set.insert(*dof_it);
                }
              }
            }
          }
        }
      }
      else
      {
        // no, this case is simpler then

        // loop over the cells in the mesh
        for (dolfin::CellIterator cell(*mesh); !cell.end(); ++cell)
        {
          // get the cell id from the mesh function
          int cellid = (*cellidmeshfunction)[(*cell).index()];
          // for each region_id we specified
          for (std::vector<int>::const_iterator id = region_ids.begin(); id != region_ids.end(); id++)
          {
            // check if this cell should be included
            if(cellid==*id)
            {
              // if yes, then get the dofs from the dofmap
              std::vector<uint> dof_vec = (*(*(*system_).functionspace())[fieldindex]).dofmap().cell_dofs((*cell).index());
              for (std::vector<uint>::const_iterator dof_it = dof_vec.begin(); dof_it < dof_vec.end(); dof_it++)
              {
                // and insert them into the set of dof for this field
                dof_set.insert(*dof_it);
              }
            }
          }
        }

      }
      // we now have a set of dofs for this field over all regions and components requested...
      // put this into the vector of dofs for all fields (field indices, unlike cell indices, should be unique so no need for a set anymore)
      for (boost::unordered_set<uint>::const_iterator dof_it = dof_set.begin(); dof_it != dof_set.end(); dof_it++)
      {
        // while we do this, check that the dofs are in the ownership range (i.e. filter out ghost nodes)
        if ((*dof_it >= ownership_range.first) && (*dof_it < ownership_range.second))
        {
          child_indices.push_back(*dof_it);
        }
      }
    }
    else
    {
      // no region ids
      // but we might still have components, so check
      buffer.str(""); buffer << optionpath << "/field[" << i << "]/components";
      if (Spud::have_option(buffer.str()))
      {
        // yes, get the components
        std::vector< int > components;
        serr = Spud::get_option(buffer.str(), components); spud_err(buffer.str(), serr);
        
        // check the user hasn't selected a component that doesn't exist
        std::vector<int>::iterator max_comp_it = std::max_element(components.begin(), components.end());
        assert(*max_comp_it < (*(*(*system_).functionspace())[fieldindex]).element().num_sub_elements());

        // loop over the components
        for (std::vector<int>::const_iterator comp = components.begin(); comp != components.end(); comp++)
        {
          // grabbing the dofs from the dofmap
          boost::unordered_set<uint> dof_set = (*(*(*(*system_).functionspace())[fieldindex])[*comp]).dofmap().dofs();
          for (boost::unordered_set<uint>::const_iterator dof_it = dof_set.begin(); dof_it != dof_set.end(); dof_it++)
          {
            // and inserting them into the vector (component dof indices should be unique, so again no need for a set)
            // check first that we own these nodes
            if ((*dof_it >= ownership_range.first) && (*dof_it < ownership_range.second))
            {
              child_indices.push_back(*dof_it);
            }
          }
        }
      }
      else
      {
        // no region ids or components specified... the simplest case
        // just grab the dofs from the dofmap
        boost::unordered_set<uint> dof_set = (*(*(*system_).functionspace())[fieldindex]).dofmap().dofs();
        for (boost::unordered_set<uint>::const_iterator dof_it = dof_set.begin(); dof_it != dof_set.end(); dof_it++)
        {
          // and insert them into the vector if we own them
          if ((*dof_it >= ownership_range.first) && (*dof_it < ownership_range.second))
          {
            child_indices.push_back(*dof_it);
          }
        }
      }
    }
  }
  // phew, that's over!
  // FIXME: tidy up the above mess!

  // sort the vector of indices over all fields (not inherited indices!)
  std::sort(child_indices.begin(), child_indices.end());

  // turn it into a simpler structure for the IS allocation
  PetscInt n=child_indices.size();
  PetscInt *indices;

  PetscMalloc(n*sizeof(PetscInt), &indices);
 
  uint ind = 0;
  if(parent_indices)
  {
    // we have been passed a list of parent indices... our child indices must be
    // a subset of this list and indexed into it so let's do that now...
    uint p_size = (*parent_indices).size();
    uint p_ind = 0;
    for (std::vector<uint>::const_iterator c_it = child_indices.begin(); c_it != child_indices.end(); c_it++)
    {
      while ((*parent_indices)[p_ind] != *c_it)
      {
        p_ind++;
        if (p_ind == p_size)
        {
          dolfin::error("Fieldsplit is not a subset of a parent fieldsplit.");
        }
      }
      indices[ind] = p_ind;
      ind++; 
      p_ind++;  // indices shouldn't be repeated so increment
    } 
    // these will probably not be equal
    assert(ind<=n);
    n = ind;
  }
  else
  {
    for (std::vector<uint>::const_iterator ind_it = child_indices.begin(); ind_it != child_indices.end(); ind_it++)
    {
      indices[ind] = *ind_it;
      ind++;
    }
    // these should be equal
    assert(ind==n);
  }

  // create the index set
  IS is;
  perr = ISCreateGeneral(PETSC_COMM_WORLD, n, indices, &is); CHKERRV(perr);
  perr = ISView(is, PETSC_VIEWER_STDOUT_SELF); CHKERRV(perr);
  perr = PCFieldSplitSetIS(pc, is); CHKERRV(perr);

  // free up the indices
  PetscFree(indices);
   
}

void SpudSolverBucket::base_fill_()
{
  // A buffer to put option paths in
  std::stringstream buffer;
  Spud::OptionError serr;

  // What is this field called?
  buffer.str(""); buffer << optionpath() << "/name";
  serr = Spud::get_option(buffer.str(), name_); spud_err(buffer.str(), serr);

  // Get the field type
  buffer.str(""); buffer << optionpath() << "/type/name";
  serr = Spud::get_option(buffer.str(), type_); spud_err(buffer.str(), serr);

  // Get the tolerances
  buffer.str(""); buffer << optionpath() << "/type/relative_error";
  serr = Spud::get_option(buffer.str(), rtol_); spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/absolute_error";
  serr = Spud::get_option(buffer.str(), atol_, 1.e-50); spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/solution_error";
  serr = Spud::get_option(buffer.str(), stol_, 1.e-8); spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/max_iterations";
  serr = Spud::get_option(buffer.str(), maxits_); spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/min_iterations";
  serr = Spud::get_option(buffer.str(), minits_, 0); spud_err(buffer.str(), serr);

  buffer.str(""); buffer << optionpath() << "/type/max_function_evaluations";
  serr = Spud::get_option(buffer.str(), maxfes_, 10000); spud_err(buffer.str(), serr);

}

void SpudSolverBucket::forms_fill_()
{
  // A buffer to put option paths (and strings) in
  std::stringstream buffer;
  Spud::OptionError serr;
   
  buffer.str(""); buffer << optionpath() << "/type/form";
  int nforms = Spud::option_count(buffer.str());
  for (uint i = 0; i < nforms; i++)
  {
    buffer.str(""); buffer << optionpath() << "/type/form[" << i << "]";
    std::string formoptionpath = buffer.str();

    std::string formname;
    buffer.str(""); buffer << formoptionpath << "/name";
    serr = Spud::get_option(buffer.str(), formname); spud_err(buffer.str(), serr);

    Form_ptr form = ufc_fetch_form((*system_).name(), name(), type(), formname, (*system_).functionspace());
    register_form(form, formname, formoptionpath);

    // Loop over the functions requested by the form and hook up (potentially currently null) pointers
    uint ncoeff = (*form).num_coefficients();
    for (uint i = 0; i < ncoeff; i++)
    {
      std::string uflsymbol = (*form).coefficient_name(i);
      GenericFunction_ptr function = (*system_).grab_uflsymbol(uflsymbol);
      if (!function)
      {
        // the function isn't initialised so either we haven't reached it yet or
        // we got to it but weren't able to initialise it because we didn't have its
        // function space available yet... let's see if it's a function
        std::string functionname = (*system_).fetch_uflname(uflsymbol);
        // this only checks the coefficient option path because fields get their functionspace
        // as a subfunctionspace from the system (mixed?) functionspace
        if (Spud::have_option((*dynamic_cast<SpudSystemBucket*>(system_)).optionpath()+"/coefficient::"+functionname+"/type::Function"))
        {
          // yes, it's a coefficient function... so let's take this opportunity to register
          // its functionspace
          if (!(*system_).contains_coefficientspace(functionname))
          {
            FunctionSpace_ptr coefficientspace;
            coefficientspace = ufc_fetch_coefficientspace((*system_).name(), name(), functionname, (*system_).mesh());
            (*system_).register_coefficientspace(coefficientspace, functionname);
          }
        }
      }

    }

  }

  if (type()=="SNES") 
  {
    linear_      = fetch_form("Residual");
    bilinear_    = fetch_form("Jacobian");
    if (contains_form("JacobianPC"))
    {
      bilinearpc_ = fetch_form("JacobianPC");
    }
    // otherwise bilinearpc_ is null
    // residual_ is a null pointer for SNES

  }
  else if (type()=="Picard")
  {
    linear_      = fetch_form("Linear");
    bilinear_    = fetch_form("Bilinear");
    if (contains_form("BilinearPC"))
    {
      bilinearpc_ = fetch_form("BilinearPC");
    }
    // otherwise bilinearpc_ is null
    residual_   = fetch_form("Residual");

  }
  else
  {
    dolfin::error("Unknown solver type.");
  }
     
}

// Register a functional in the function
void SpudSolverBucket::register_form(Form_ptr form, std::string name, std::string optionpath)
{
  // First check if a form with this name already exists
  Form_it f_it = forms_.find(name);
  if (f_it != forms_.end())
  {
    // if it does, issue an error
    dolfin::error("Form named \"%s\" already exists in function.", name.c_str());
  }
  else
  {
    // if not then insert it into the maps
    forms_[name]            = form;
    form_optionpaths_[name] = optionpath;
  }
}

// Return a string describing the contents of the spudfunction
const std::string SpudSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Describe the contents of the form_optionpaths_ map
const std::string SpudSolverBucket::forms_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = form_optionpaths_.begin(); s_it != form_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

