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
  int numCells() const { return numCells_; }

  double totalTDCWidth() const;
  double maxTDCWidth() const;

  const TDCHitWPPtrCollection& hits() const { return hits_; }

  explicit WireCluster(const TDCHitWPPtrCollection& hits);

private:
  TDCHitWPPtrCollection hits_;
  int plane_;
  double centralCell_;
  int numCells_;
};

// clusters in one plane
typedef std::vector<WireCluster> WireClusterCollection;

// clusters[plane]
typedef std::vector<WireClusterCollection> ClustersByPlane;

// The clusterization, for PC and DC hits
ClustersByPlane constructPlaneClusters(int maxPlaneNumber, const TDCHitWPPtrCollection& hits);

// Combine PCs and DCs together
ClustersByPlane globalPlaneClusters(const ClustersByPlane& pcClusters, const ClustersByPlane& dcClusters);

TDCHitWPPtr maxTDCWidthHit(const WireCluster& cluster);
// the caller must protect against emtpy planeClusters
TDCHitWPPtr maxTDCWidthHit(const WireClusterCollection& planeClusters);

std::ostream& operator<<(std::ostream& os, const WireCluster& cl);
std::ostream& operator<<(std::ostream& os, const WireClusterCollection& coll);
std::ostream& operator<<(std::ostream& os, const ClustersByPlane& cbp);

#endif/*WireCluster_h*/
