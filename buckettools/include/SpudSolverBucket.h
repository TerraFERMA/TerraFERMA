
#ifndef __SPUD_SOLVERBUCKET_H
#define __SPUD_SOLVERBUCKET_H

#include "BoostTypes.h"
#include "SolverBucket.h"
#include <dolfin.h>

namespace buckettools
{
  
  // The SpudSolverBucket class is a derived class of the solver that populates the
  // data structures within a solver using the spud options parser (and assumes
  // the structure of the buckettools schema)
  class SpudSolverBucket : public SolverBucket
  {
  // only accessible by this class
  private:

    // supplement the base class with an optionpath
    std::string optionpath_;

    // supplementary to the solver base data store the optionpaths for the forms_
    std::map< std::string, std::string > form_optionpaths_;
    
    // fill in the base data
    void base_fill_();
 
    // fill in the information related to the forms of this solver
    void forms_fill_();

    void ksp_fill_(const std::string &optionpath, KSP &ksp)
    { ksp_fill_(optionpath, ksp, NULL); }

    void ksp_fill_(const std::string &optionpath, KSP &ksp, const std::vector<uint>* parent_indices);

    void pc_fieldsplit_fill_(const std::string &optionpath, PC &pc, const std::vector<uint>* parent_indices);

    void pc_fieldsplit_by_field_fill_(const std::string &optionpath, PC &pc, 
                                      std::vector<uint> &child_indices,
                                      const std::vector<uint>* parent_indices);

  public:
    
    // Specific constructor
    SpudSolverBucket(std::string optionpath, SystemBucket* system);
    
    // Default destructor (virtual so it calls base class as well)
    virtual ~SpudSolverBucket();

    // Fill the solver assuming the buckettools schema
    void fill();

    // Register a form in the solver (with a spud optionpath)
    void register_form(Form_ptr form, std::string name, std::string optionpath);

    // Return a string object describing the solver
    const std::string str() const
    { return str(0); }

    // Return a string object describing the solver
    const std::string str(int indent) const;

    // Return a string object describing the forms
    const std::string forms_str() const
    { return forms_str(0); }

    // Return a string object describing the forms
    const std::string forms_str(int indent) const;

    // Return the base optionpath for this function
    const std::string optionpath() const
    { return optionpath_; }
    
  };
 
  // Define a boost shared pointer for a spud function
  typedef boost::shared_ptr< SpudSolverBucket > SpudSolverBucket_ptr;

}
#endif
