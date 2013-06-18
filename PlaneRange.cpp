// Andrei Gaponenko, 2013

#include "PlaneRange.h"

PlaneRange findPlaneRange(const ClustersByPlane& cp) {
  PlaneRange::Segments res;
  unsigned currentPlane = 0;

  while(++currentPlane < cp.size()) {

    if(!cp[currentPlane].empty()) { // start a new segment

      const int pstart = currentPlane;
      while(++currentPlane < cp.size() && !cp[currentPlane].empty())
        {}

      const int pend = currentPlane - 1;
      res.emplace_back(pstart, pend);
    }
  }

  return PlaneRange(res);
}
