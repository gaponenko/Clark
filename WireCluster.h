// Andrei Gaponenko, 2013

#ifndef WireCluster_h
#define WireCluster_h

#include <vector>
#include <map>
#include <ostream>

#include "TDCHitWP.h"

class WireCluster {
public:

  int plane() const { return plane_; }
  double centralCell() const { return centralCell_; }
  double numCells() const { return numCells_; }

  const TDCHitWPPtrCollection& hits() const { return hits_; }

  explicit WireCluster(const TDCHitWPPtrCollection& hits);

private:
  TDCHitWPPtrCollection hits_;
  int plane_;
  double centralCell_;
  int numCells_;
};

// plane number => clusters
typedef std::vector<WireCluster> WireClusterCollection;
typedef std::map<int, WireClusterCollection> ClustersByPlane;

std::ostream& operator<<(std::ostream& os, const WireCluster& cl);
std::ostream& operator<<(std::ostream& os, const WireClusterCollection& coll);
std::ostream& operator<<(std::ostream& os, const ClustersByPlane& cbp);

#endif/*WireCluster_h*/
