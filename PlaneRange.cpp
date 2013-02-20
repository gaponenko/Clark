// Andrei Gaponenko, 2013

#include "PlaneRange.h"

PlaneRange findPlaneRange(const ClustersByPlane& cp) {
  PlaneRange res;
  typedef ClustersByPlane::const_iterator Iter;

  int numHitPlanes(0);
  for(unsigned  i = 0; i<cp.size(); ++i) {
    if(!cp[i].empty()) {
      ++numHitPlanes;
      if(res.min == -1) {
        res.min = i;
        res.max = i;
      }
      else {
        if(res.min > i) {
          res.min = i;
        }
        if(i > res.max) {
          res.max= i;
        }
      }
    }
  }

  res.noGaps = ((res.max - res.min + 1) == numHitPlanes);

  return res;
}
