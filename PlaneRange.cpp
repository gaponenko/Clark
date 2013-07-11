// Andrei Gaponenko, 2013

#include "PlaneRange.h"
#include <cassert>

//================================================================
PlaneRange findRestrictedPlaneRange(const ClustersByPlane& cp, int begin, int end) {
  assert(begin >= 0);
  assert(end <= cp.size());

  PlaneRange::Segments res;
  int currentPlane = begin - 1;

  while(++currentPlane < end) {

    if(!cp[currentPlane].empty()) { // start a new segment

      const int pstart = currentPlane;
      while(++currentPlane < end && !cp[currentPlane].empty())
        {}

      const int pend = currentPlane - 1;
      res.push_back(PlaneRangeSegment(pstart, pend));
    }
  }

  return PlaneRange(res);
}

//================================================================
PlaneRange findPlaneRange(const ClustersByPlane& cp) {
  return findRestrictedPlaneRange(cp, 0, cp.size());
}

//================================================================
PlaneRange findUpstreamPlaneRange(const ClustersByPlane& cp) {
  assert(cp.size() == 57);
  return findRestrictedPlaneRange(cp, 1, 29);
}

//================================================================
PlaneRange findDownstreamPlaneRange(const ClustersByPlane& cp) {
  assert(cp.size() == 57);
  return findRestrictedPlaneRange(cp, 29, cp.size());
}

//================================================================
std::ostream& operator<<(std::ostream& os, const PlaneRange& r) {
  os<<"PlaneRange(";
  for(unsigned i=0; i<r.segments().size(); ++i) {
    os<<"["<<r.segments()[i].min<<", "<<r.segments()[i].max<<"], ";
  }
  return os<<") ";
}

//================================================================
