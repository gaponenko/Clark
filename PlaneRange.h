// Andrei Gaponenko, 2013

#ifndef PlaneRange_h
#define PlaneRange_h

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

PlaneRange findPlaneRange(const ClustersByPlane& clusters);


#endif/*PlaneRange_h*/
