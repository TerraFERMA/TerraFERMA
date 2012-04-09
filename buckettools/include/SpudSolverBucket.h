
#ifndef __SPUD_SOLVERBUCKET_H
#define __SPUD_SOLVERBUCKET_H

#include "BoostTypes.h"
#include "SolverBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  //*****************************************************************|************************************************************//
  // SpudSolverBucket class:
  //
  // The SpudSolverBucket class is a derived class of the solver that populates the
  // data structures within a solver using the spud options parser (and assumes
  // the structure of the buckettools schema)
  //*****************************************************************|************************************************************//
  class SpudSolverBucket : public SolverBucket
  {

  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//

  public:
    
    //***************************************************************|***********************************************************//
    // Constructors and destructors
    //***************************************************************|***********************************************************//

    SpudSolverBucket(const std::string &optionpath, SystemBucket* system);  // specific constructor (taking in optionpath and parent system)
    
    ~SpudSolverBucket();                                             // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill in the data in the base solver bucket

    void copy_diagnostics(SolverBucket_ptr &solver, 
                                  SystemBucket_ptr &system) const;   // copy the data necessary for the diagnostics data file(s)

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a constant string containing the optionpath for the solver bucket
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Form data access
    //***************************************************************|***********************************************************//

    void register_form(Form_ptr form, const std::string &name,       // register a form with the given name and optionpath in the solver
                                   const std::string &optionpath);   // bucket

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str() const                                    // return a string describing the contents of the solver bucket
    { return str(0); }

    const std::string str(int indent) const;                         // return an indented string describing the contents of the solver bucket

    const std::string forms_str() const                              // return a string describing the forms in the solver bucket
    { return forms_str(0); }

    const std::string forms_str(const int &indent) const;            // return an indented string describing the forms in the solver bucket

  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//

  private:                                                           // only accessible by this class

    //***************************************************************|***********************************************************//
    // Base data
    //***************************************************************|***********************************************************//

    std::string optionpath_;                                         // the optionpath for the solver bucket

    //***************************************************************|***********************************************************//
    // Pointers data
    //***************************************************************|***********************************************************//

    std::map< std::string, std::string > form_optionpaths_;          // a map from form names to form optionpaths
    
    //***************************************************************|***********************************************************//
    // Filling data (continued)
    //***************************************************************|***********************************************************//

    void fill_base_();                                               // fill the base data of the solver bucket
 
    void fill_forms_();                                              // fill the form data of the solver bucket

    void fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                  const std::string prefix)          // fill the information about a parent ksp
    { fill_ksp_(optionpath, ksp, prefix, NULL, NULL); }

    void fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                  const std::string prefix,
                                  SNES *snes)                  // fill the information about a parent ksp
    { fill_ksp_(optionpath, ksp, prefix, NULL, snes); }

    void fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                  const std::string prefix,          // fill the information about a child ksp
                                  const std::vector<uint>* parent_indices)
    { fill_ksp_(optionpath, ksp, prefix, parent_indices, NULL); }

    void fill_ksp_(const std::string &optionpath, KSP &ksp,          // fill the information about a child ksp
                         const std::string prefix, 
                         const std::vector<uint>* parent_indices,
                         SNES *snes);

    void fill_pc_fieldsplit_(const std::string &optionpath, PC &pc,  // fill the information about a fieldsplit pc
                         const std::string prefix, 
                         const std::vector<uint>* parent_indices);

    void fill_is_by_field_(const std::string &optionpath, IS &is,    // set up a petsc index set
                           std::vector<uint> &child_indices,
                           const std::vector<uint>* parent_indices,
                           const std::vector<uint>* sibling_indices);

    boost::unordered_set<uint> region_dof_set_(const std::string &optionpath, // return a vector of dofs in a particular region
                                               const boost::shared_ptr<const dolfin::GenericDofMap> dofmap);

    boost::unordered_set<uint> boundary_dof_set_(const std::string &optionpath, // return a vector of dofs in a particular region
                                                 const boost::shared_ptr<const dolfin::GenericDofMap> dofmap);

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                           // empty the derived and base class data structures

    
  };
 
  typedef boost::shared_ptr< SpudSolverBucket > SpudSolverBucket_ptr;// define a (boost shared) pointer to a spud solver bucket class

}
#endif
