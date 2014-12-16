// The difference from HistTDCWidth is that here we
// histogram only selected hits from a given time window:
// if there is more than one cluster in a plane, those
// hits are not used.  Accepted hits are histogrammed
// in bins of cluster size ("By Cluster Size").
//
// Andrei Gaponenko, 2014

#ifndef HistTDCBCSWidth_h
#define HistTDCBCSWidth_h

#include <vector>

#include "WireCluster.h"
#include "HistTDCWidth.h"
#include "DetectorGeo.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistTDCBCSWidth {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt, const ClustersByPlane& globalClusters);

  void fill(const EventClass& evt, WirePlane::DetType pt, const WireCluster& c);

private :
  static const int maxClusterSize = 4;

  const DetectorGeo *geom_;

  TH2 *clusterSizePC_;
  TH2 *clusterSizeDC_;

  std::vector<HistTDCWidth> clusterSizeBinsPC_;
  std::vector<HistTDCWidth> clusterSizeBinsDC_;

  HistTDCWidth perWireHitDistsPC_;
  HistTDCWidth perWireHitDistsDC_;
};

#endif/*HistTDCBCSWidth_h*/
