// Andrei Gaponenko, 2013

#ifndef PlaneRange_h
#define PlaneRange_h

#include "WireCluster.h"

struct PlaneRange {
  int min;
  int max;
  bool noGaps;
  PlaneRange() : min(-1), max(-1), noGaps(false) {}
};

PlaneRange findPlaneRange(const ClustersByPlane& clusters);


#endif/*PlaneRange_h*/
