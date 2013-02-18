// Andrei Gaponenko, 2013

#ifndef WireCluster_h
#define WireCluster_h

#include <vector>
#include <map>

#include "TDCHitWP.h"

struct WireCluster {
  int plane;
  std::vector<TDCHitWPPtr> hits;
  WireCluster() : plane(-1) {}
};

// plane number => clusters
typedef std::vector<WireCluster> WireClusterCollection;
typedef std::map<int, WireClusterCollection> ClustersByPlane;

#endif/*WireCluster_h*/
