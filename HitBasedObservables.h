// Andrei Gaponenko, 2014

#ifndef HitBasedObservables_h
#define HitBasedObservables_h

#include <vector>
#include "WireCluster.h"

//================================================================
class HitBasedObservables {
public:
  explicit HitBasedObservables(const ClustersByPlane& protonGlobalClusters);
  unsigned dnCPlanes() const { return dnCPlanes_; }
  unsigned dnCWires() const { return dnCWires_; }

  typedef std::vector<int> ClusterSizeList;
  ClusterSizeList clusterSize() const { return clusterSize_; }

private :
  unsigned dnCPlanes_; // num contiguous planes hit. Zero if no hits in PC7
  unsigned dnCWires_; // sum of largest cluster sizes per contiguous plane
  ClusterSizeList clusterSize_;
};

#endif/*HitBasedObservables_h*/
