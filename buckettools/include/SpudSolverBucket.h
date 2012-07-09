
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

    void initialize();                                               // initialize the solvers and tensors

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

    void fill_subforms_(const std::string &optionpath, 
                        const std::string &prefix="");               // fill the form data of a section of the solver bucket

    void fill_solverforms_(const std::string &optionpath, 
                           const std::string &prefix="");            // fill the form data of a linear solver

    void fill_ksp_(const std::string &optionpath, KSP &ksp, 
                                  const std::string prefix)          // fill the information about a parent ksp
    { fill_ksp_(optionpath, ksp, prefix, NULL); }

    void fill_ksp_(const std::string &optionpath, KSP &ksp,          // fill the information about a child ksp
                         const std::string prefix, 
                         const std::vector<uint>* parent_indices);

    void fill_pc_(const std::string &optionpath, PC &pc,             // fill the information about a pc
                         const std::string prefix, 
                         const std::vector<uint>* parent_indices);

    void fill_pc_fieldsplit_(const std::string &optionpath, PC &pc,  // fill the information about a fieldsplit pc
                         const std::string prefix, 
                         const std::vector<uint>* parent_indices);

    void fill_is_by_field_(const std::string &optionpath, IS &is,    // set up a petsc index set
                           std::vector<uint> &child_indices,
                           const std::vector<uint>* parent_indices,
                           const std::vector<uint>* sibling_indices);

    boost::unordered_set<uint> field_dof_set_(const std::string &optionpath,
                                              const FunctionSpace_ptr functionspace,
                                              const std::vector<int>* components,
                                              const std::vector<int>* region_ids,
                                              const std::vector<int>* boundary_ids,
                                              const uint parent_component=0,
                                              uint rank=0);          // set up a dof set based on a field

    boost::unordered_set<uint> cell_dof_set_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                             const std::vector<int>* region_ids);
                                                                     // set up a dof set over cells

    boost::unordered_set<uint> facet_dof_set_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                              const std::vector<int>* boundary_ids);
                                                                     // set up a dof set over facets

    void fill_nullspace_(const std::string &optionpath, 
                         MatNullSpace &SP,
                         const std::vector<uint>* parent_indices);   // fill a petsc null space object

    void fill_constraints_();                                        // fill constraints on snes vi

    void fill_bound_(const std::string &optionpath, 
                         PETScVector_ptr &bound, const double &background_value);   // fill petsc vectors containing the bounds for snes vi

    void fill_values_by_field_(const std::string &optionpath,    // fill a vector describing e.g. a null space
                                   PETScVector_ptr values,
                                   const double &value,
                                   const std::vector<uint>* parent_indices,
                                   const std::vector<uint>* sibling_indices);

    boost::unordered_map<uint, double> field_value_map_(const std::string &optionpath,
                                                     const FunctionSpace_ptr functionspace,
                                                     const std::vector<int>* components,
                                                     const std::vector<int>* region_ids,
                                                     const std::vector<int>* boundary_ids,
                                                     Expression_ptr value_exp, const double *value_const,
                                                     const uint parent_component=0,
                                                     uint rank=0,    // set up a map describing e.g. the null space values for a field
                                                     uint exp_index=0);

    boost::unordered_map<uint, double> cell_value_map_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                    const std::vector<int>* region_ids,
                                                    Expression_ptr value_exp, const double *value_const,
                                                    const uint &exp_index);// set up a map describing e.g. the null space values over cells

    boost::unordered_map<uint, double> facet_value_map_(const boost::shared_ptr<const dolfin::GenericDofMap> dofmap,
                                                     const std::vector<int>* boundary_ids,
                                                     Expression_ptr value_exp, const double *value_const,
                                                     const uint &exp_index);// set up a map describing e.g. the  null space values over a boundary

    void is_by_field_restrictions_(const std::string &optionpath,
                                   std::vector<int>* &components,
                                   std::vector<int>* &region_ids,
                                   std::vector<int>* &boundary_ids,
                                   const std::string &fieldrank,
                                   const int &fieldsize);            // set up the restrictions on an IS by field

    void destroy_is_field_restrictions_(std::vector<int>* &components,// destroy the objects describing any restrictions on an IS by field
                                        std::vector<int>* &region_ids,
                                        std::vector<int>* &boundary_ids);

    void restrict_is_indices_(std::vector<uint> &indices,            // restrict the indices describing an IS based on the parent,
                              const std::vector<uint>* parent_indices,// sibling and parallel ownership
                              const std::vector<uint>* sibling_indices);

    void initialize_tensors_();                                      // fill the tensor data structures of the solver bucket

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                           // empty the derived and base class data structures

    
  };
 
  typedef boost::shared_ptr< SpudSolverBucket > SpudSolverBucket_ptr;// define a (boost shared) pointer to a spud solver bucket class

}
#endif
