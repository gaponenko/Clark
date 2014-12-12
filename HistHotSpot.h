// Various plots to understand the "hot spot" in the hit multiplicity
// distribution.
//
// Andrei Gaponenko, 2014

#ifndef HistHotSpot_h
#define HistHotSpot_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistHotSpot {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt, const ClustersByPlane& globalPlaneClusters);

private :
  // One entry per event
  TH1* hitTimeWidest_;

  TH2* tdcWidthWidestHit8vs7_;
  TH2* posWidestHit8vs7_;

  TH2* clusterSize8vs7_;

  TH1* numUnassignedDCHits_;

  // Multiple entries per event
  TH1 *hitTimeAll78_;
  std::vector<TH1*> hitTimeAllPlane_; //PC7, PC8

  TH1 *tdcWidthAll78_;
  std::vector<TH1*> tdcWidthAllPlane_;

  std::vector<TH1*> dt_; //PC7, PC8

  // fraction of total in N widest hits vs N? (Or vs hit fraction: both axes use [0:1])
  /// can we find a translation-invariant measure?
  // Kyle: plot (xn1 - xn)
  // mean vs median?  visualize "one broad, many narrow".

  // correlation with pc5/6? Other upstream?
};

#endif/*HistHotSpot_h*/
