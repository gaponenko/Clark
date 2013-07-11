// Andrei Gaponenko, 2013

#ifndef PlaneRange_h
#define PlaneRange_h

#include <ostream>

#include "WireCluster.h"

struct PlaneRangeSegment {
  int min;
  int max;
  PlaneRangeSegment(int n1, int n2) : min(n1), max(n2) {}
};

class PlaneRange {
public:
  typedef std::vector<PlaneRangeSegment> Segments;
  const Segments& segments() const { return segments_; }
  int min() const { return segments_.at(0).min; }
  int max() const { return segments_.at(segments_.size()-1).max; }
  bool noGaps() const { return segments_.size() == 1; }

  PlaneRange(const Segments& s) : segments_(s) {}
private:
  std::vector<PlaneRangeSegment> segments_;
};

// end is one-past-end, a la STL
PlaneRange findRestrictedPlaneRange(const ClustersByPlane& clusters, int begin, int end);

PlaneRange findPlaneRange(const ClustersByPlane& clusters);
PlaneRange findUpstreamPlaneRange(const ClustersByPlane& clusters);
PlaneRange findDownstreamPlaneRange(const ClustersByPlane& clusters);


std::ostream& operator<<(std::ostream& os, const PlaneRange& r);

#endif/*PlaneRange_h*/
