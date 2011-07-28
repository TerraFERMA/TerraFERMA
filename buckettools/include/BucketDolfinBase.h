
#ifndef __BUCKETDOLFIN_BASE_H
#define __BUCKETDOLFIN_BASE_H

#include <dolfin.h>

namespace buckettools
{

  class Side : public dolfin::SubDomain
  {
  private:
    uint component_;
    double side_;

  public:
    Side(const uint &component, const double &side);

    ~Side();

    bool inside(const dolfin::Array<double>& x, bool on_boundary) const; 

  };

}

#endif

