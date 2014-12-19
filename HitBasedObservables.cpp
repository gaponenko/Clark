// Andrei Gaponenko, 2014

#include "HitBasedObservables.h"

template<class ClusterCmp>
HitBasedObservables<ClusterCmp>::HitBasedObservables(const ClustersByPlane& protonGlobalClusters)
  : dnCPlanes_(0)
  , dnCWires_(0)
{
  // Start at PC7
  for(int iplane = 29; (iplane < protonGlobalClusters.size()) && !protonGlobalClusters.at(iplane).empty(); ++iplane) {
    const WireClusterCollection& clusters = protonGlobalClusters[iplane];

    // Find the best cluster
    WireClusterCollection::const_iterator best = clusters.begin();
    for(WireClusterCollection::const_iterator current = ++clusters.begin(); current != clusters.end(); ++current) {
      if(cmp_(*best, *current)) {
        best = current;
      }
    }
    clusters_.push_back(*best);
  }

  dnCPlanes_ = clusters_.size();

  dnCWires_ = 0;
  for(int i=0; i<clusters_.size(); ++i) {
    dnCWires_ += clusters_[i].numCells();
  }
}

template HitBasedObservables<WireClusterCmpNumCells>::HitBasedObservables(const ClustersByPlane&);
template HitBasedObservables<WireClusterCmpTDCWidth>::HitBasedObservables(const ClustersByPlane&);
