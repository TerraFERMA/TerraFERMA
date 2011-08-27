
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

    SpudSolverBucket(std::string optionpath, SystemBucket* system);  // specific constructor (taking in optionpath and parent system)
    
    ~SpudSolverBucket();                                             // default destructor

    //***************************************************************|***********************************************************//
    // Filling data
    //***************************************************************|***********************************************************//

    void fill();                                                     // fill in the data in the base solver bucket

    //***************************************************************|***********************************************************//
    // Base data access
    //***************************************************************|***********************************************************//

    const std::string optionpath() const                             // return a constant string containing the optionpath for the solver bucket
    { return optionpath_; }

    //***************************************************************|***********************************************************//
    // Form data access
    //***************************************************************|***********************************************************//

    void register_form(Form_ptr form, std::string name,              // register a form with the given name and optionpath in the solver
                                          std::string optionpath);   // bucket

    //***************************************************************|***********************************************************//
    // Output functions
    //***************************************************************|***********************************************************//

    const std::string str() const                                    // return a string describing the contents of the solver bucket
    { return str(0); }

    const std::string str(int indent) const;                         // return an indented string describing the contents of the solver bucket

    const std::string forms_str() const                              // return a string describing the forms in the solver bucket
    { return forms_str(0); }

    const std::string forms_str(int indent) const;                   // return an indented string describing the forms in the solver bucket

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

    void base_fill_();                                               // fill the base data of the solver bucket
 
    void forms_fill_();                                              // fill the form data of the solver bucket

    void ksp_fill_(const std::string &optionpath, KSP &ksp)          // fill the information about a parent ksp
    { ksp_fill_(optionpath, ksp, NULL); }

    void ksp_fill_(const std::string &optionpath, KSP &ksp,          // fill the information about a child ksp
                         const std::vector<uint>* parent_indices);

    void pc_fieldsplit_fill_(const std::string &optionpath, PC &pc,  // fill the information about a fieldsplit pc
                         const std::vector<uint>* parent_indices);

    void is_by_field_fill_(const std::string &optionpath, IS &is,    // set up a petsc index set
                           std::vector<uint> &child_indices,
                           const std::vector<uint>* parent_indices);

    //***************************************************************|***********************************************************//
    // Emptying data
    //***************************************************************|***********************************************************//

    void empty_();                                           // empty the derived and base class data structures

    
  };
 
  typedef boost::shared_ptr< SpudSolverBucket > SpudSolverBucket_ptr;// define a (boost shared) pointer to a spud solver bucket class

}
#endif
