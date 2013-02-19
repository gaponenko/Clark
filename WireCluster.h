// Andrei Gaponenko, 2013

#ifndef WireCluster_h
#define WireCluster_h

#include <vector>
#include <map>
#include <ostream>

#include "TDCHitWP.h"

struct WireCluster {
  int plane;
  std::vector<TDCHitWPPtr> hits;
  WireCluster() : plane(-1) {}
};

// plane number => clusters
typedef std::vector<WireCluster> WireClusterCollection;
typedef std::map<int, WireClusterCollection> ClustersByPlane;

std::ostream& operator<<(std::ostream& os, const WireCluster& cl);
std::ostream& operator<<(std::ostream& os, const WireClusterCollection& coll);
std::ostream& operator<<(std::ostream& os, const ClustersByPlane& cbp);

#endif/*WireCluster_h*/
