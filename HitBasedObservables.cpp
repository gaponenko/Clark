// Andrei Gaponenko, 2014

#include "HitBasedObservables.h"

HitBasedObservables::HitBasedObservables(const ClustersByPlane& protonGlobalClusters)
  : dnCPlanes_(0)
  , dnCWires_(0)
{
  // Require a hit in PC7
  if(!protonGlobalClusters.at(29).empty()) {

    dnCPlanes_ = 28;

    do {
      ++dnCPlanes_;
      int maxclustersize = 0;
      for(int i=0; i<protonGlobalClusters.at(dnCPlanes_).size(); ++i) {
        maxclustersize = std::max(maxclustersize, protonGlobalClusters[dnCPlanes_][i].numCells());
      }
      dnCWires_ += maxclustersize;
      clusterSize_.push_back(maxclustersize);
    } while(!protonGlobalClusters.at(1+dnCPlanes_).empty());

    dnCPlanes_ -= 28;
  }
}
