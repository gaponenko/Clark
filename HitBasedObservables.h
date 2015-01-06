// Andrei Gaponenko, 2014

#ifndef HitBasedObservables_h
#define HitBasedObservables_h

#include <vector>
#include <string>
#include "WireCluster.h"

class TH2;
class HistogramFactory;
class DetectorGeo;
class ConfigFile;

//================================================================
// Histograms that can optionally be filled during the
// HitBasedObservables computation.
class HistHitBasedAmbiguities {
public:
  TH2 *hClusterMultiplicity_;
  TH2 *hDiffNWires_;

  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);
};

//================================================================
template<class ClusterCmp>
class HitBasedObservables {
public:
  explicit HitBasedObservables(const ClustersByPlane& protonGlobalClusters, HistHitBasedAmbiguities *hh=0);

  // count of contiguous planes hit
  unsigned dnCPlanes() const { return dnCPlanes_; }
  // total sum of "contiguous" wires (cluster sizes)
  unsigned dnCWires() const { return dnCWires_; }

  // Clusters selected for the hit-based analysis, one per plane,
  // starting with PC7 at [0].
  WireClusterCollection clusters() const { return clusters_; }

private :
  ClusterCmp cmp_;
  unsigned dnCPlanes_; // num contiguous planes hit. Zero if no hits in PC7
  unsigned dnCWires_; // sum of largest cluster sizes per contiguous plane
  WireClusterCollection clusters_;
};

typedef HitBasedObservables<WireClusterCmpNumCells> HitBasedObservablesMaxSize;
typedef HitBasedObservables<WireClusterCmpTDCWidth> HitBasedObservablesMaxWidth;

#endif/*HitBasedObservables_h*/
