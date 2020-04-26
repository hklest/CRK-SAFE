#ifndef __GEOM_SURFACE_EFF_H__
#define __GEOM_SURFACE_EFF_H__

#include "Eff.h"

class Geom_Surface_Eff : public Eff
{
public:
  Geom_Surface_Eff (double EFF) {eff = EFF;}
  virtual ~Geom_Surface_Eff() {}

  // This efficiency has no wavelength dependence
  double e(double lambda)
  {
    eff = .8333; // GEMs are 83% material and 17% holes, holes have 0 efficiency.
    return eff;
  }

protected:
  double eff;
};
	
#endif /* __GEOM_SURFACE_EFF_H__ */
