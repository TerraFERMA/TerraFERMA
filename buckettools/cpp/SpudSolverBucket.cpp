
#include "BoostTypes.h"
#include "SpudSolverBucket.h"
#include <dolfin.h>
#include <string>
#include <spud>
#include "SystemsWrapper.h"
#include "SpudSystem.h"
#include "SpudBase.h"
#include "SpudBucket.h"
#include "petscsnes.h"

using namespace buckettools;

// Specific constructor
SpudSolverBucket::SpudSolverBucket(std::string optionpath, System* system) : optionpath_(optionpath), SolverBucket(system)
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
    buffer.str(""); buffer << name() << "_";
    perr = SNESSetOptionsPrefix(snes_, buffer.str().c_str()); CHKERRV(perr);

    // we always have at least one ksp
    perr = SNESGetKSP(snes_, &ksp_); CHKERRV(perr);
    
    buffer.str(""); buffer << optionpath() << "/type/linear_solver";
    ksp_fill_(buffer.str(), ksp_);
  }
  else if (type()=="Picard")
  {
;
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
    pc_fieldsplit_by_field_fill_(buffer.str(), pc, child_indices[i], parent_indices);
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

  // a vecto of global indices
  std::vector< uint > indices_vector;
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
        std::vector<int>::iterator max_comp_it = std::max(components.begin(), components.end());
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
          indices_vector.push_back(*dof_it);
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
        std::vector<int>::iterator max_comp_it = std::max(components.begin(), components.end());
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
              indices_vector.push_back(*dof_it);
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
            indices_vector.push_back(*dof_it);
          }
        }
      }
    }
  }
  // phew, that's over!
  // FIXME: tidy up the above mess!

  // sort the vector of indices over all fields (not inherited indices!)
  std::sort(indices_vector.begin(), indices_vector.end());

  // turn it into a simpler structure for the IS allocation
  PetscInt n=indices_vector.size();
  PetscInt *indices;

  PetscMalloc(n*sizeof(PetscInt), &indices);
  uint ind = 0;
  for (std::vector<uint>::const_iterator ind_it = indices_vector.begin(); ind_it != indices_vector.end(); ind_it++)
  {
    indices[ind] = *ind_it;
    ind++;
  }

  // these should be equal
  assert(ind==n);

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
      GenericFunction_ptr function = (*system_).fetch_uflsymbol(uflsymbol);
      if (!function)
      {
        // the function isn't initialised so either we haven't reached it yet or
        // we got to it but weren't able to initialise it because we didn't have its
        // function space available yet... let's see if it's a function
        std::string functionname = (*system_).fetch_uflname(uflsymbol);
        if (Spud::have_option((*dynamic_cast<SpudSystem*>(system_)).optionpath()+"/coefficient::"+functionname+"/type::Function"))
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

      (*form).set_coefficient(uflsymbol, function);
    }
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
std::string SpudSolverBucket::str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');
  s << indentation << "SolverBucket " << name() << " (" << optionpath() << ")" << std::endl;
  indent++;
  s << forms_str(indent);
  return s.str();
}

// Describe the contents of the form_optionpaths_ map
std::string SpudSolverBucket::forms_str(int indent) const
{
  std::stringstream s;
  std::string indentation (indent*2, ' ');

  for ( string_const_it s_it = form_optionpaths_.begin(); s_it != form_optionpaths_.end(); s_it++ )
  {
    s << indentation << "Form " << (*s_it).first << " (" << (*s_it).second  << ")" << std::endl;
  }

  return s.str();
}

