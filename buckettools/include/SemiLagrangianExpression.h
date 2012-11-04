
#ifndef __SEMI_LAGRANGIAN_EXPRESSION_H
#define __SEMI_LAGRANGIAN_EXPRESSION_H

#include "BoostTypes.h"
#include "Bucket.h"
#include "SystemBucket.h"
#include <dolfin.h>

namespace buckettools
{
  struct lt_point
  {
    bool operator() (const dolfin::Point& p1, const dolfin::Point& p2) const
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        if (p1[i] < (p2[i] - DOLFIN_EPS))
          return true;
        else if (p1[i] > (p2[i] + DOLFIN_EPS))
          return false;
      }

      return false;
    }
  };
  
  typedef std::map<dolfin::Point, std::vector<int>, lt_point> point_map;
  typedef point_map::iterator point_iterator;


  //*****************************************************************|************************************************************//
  // SemiLagrangianExpression class:
  //
  // The SemiLagrangianExpression class describes a derived dolfin Expression class that overloads
  // the eval function using a user defined data.
  //*****************************************************************|************************************************************//
  class SemiLagrangianExpression : public dolfin::Expression
  {
  
  //*****************************************************************|***********************************************************//
  // Publicly available functions
  //*****************************************************************|***********************************************************//
  
  public:                                                            // available to everyone
  
    SemiLagrangianExpression(const Bucket *bucket, const double_ptr time, 
                             const std::pair< std::string, std::pair< std::string, std::string > > &function,
                             const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                             const std::pair< std::string, std::pair< std::string, std::string > > &outside);

    SemiLagrangianExpression(const uint &dim,
                             const Bucket *bucket, const double_ptr time, 
                             const std::pair< std::string, std::pair< std::string, std::string > > &function,
                             const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                             const std::pair< std::string, std::pair< std::string, std::string > > &outside);

    SemiLagrangianExpression(const uint &dim0, const uint &dim1,
                             const Bucket *bucket, const double_ptr time, 
                             const std::pair< std::string, std::pair< std::string, std::string > > &function,
                             const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                             const std::pair< std::string, std::pair< std::string, std::string > > &outside);

    SemiLagrangianExpression(const std::vector<uint> &value_shape,
                             const Bucket *bucket, const double_ptr time, 
                             const std::pair< std::string, std::pair< std::string, std::string > > &function,
                             const std::pair< std::string, std::pair< std::string, std::string > > &velocity,
                             const std::pair< std::string, std::pair< std::string, std::string > > &outside);

    ~SemiLagrangianExpression();
    
    void eval(dolfin::Array<double>& values, 
              const dolfin::Array<double>& x, 
              const ufc::cell &cell) const;
    
    void init();
  
  //*****************************************************************|***********************************************************//
  // Private functions
  //*****************************************************************|***********************************************************//
  
  private:                                                           // only available to this class
    
    const Bucket* bucket() const
    { return bucket_; }
    
    const double_ptr time() const
    {
      return time_;
    }

    const Bucket *bucket_;
    
    const double_ptr time_;

    bool initialized_;

    std::pair< std::string, std::pair< std::string, std::string > > 
                                        funcname_, velname_, outname_;
    
    GenericFunction_ptr vel_, oldvel_;
    GenericFunction_ptr func_, out_;

    Mesh_ptr mesh_;
    uint dim_;

    double *xstar_;

    dolfin::CellIterator *dolfincellit_;
    dolfin::UFCCell *ufccellstar_;

    dolfin::Array<double> *v_, *oldv_, *vstar_;

    dolfin::MeshFunction< point_map > *cellcache_;
    
    const bool findpoint_(const dolfin::Array<double>& x, 
                          const ufc::cell &cell) const;

    void findvstar_(const dolfin::Array<double>& x, 
                    const ufc::cell &cell) const;

    const bool checkpoint_(const int &index, 
                           point_iterator &p_it, 
                           point_map &points, 
                           const dolfin::Point &lp) const;

  };

}

#endif

