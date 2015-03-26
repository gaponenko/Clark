// Andrei Gaponenko, 2015

#ifndef HistMuCapMuonRange_h
#define HistMuCapMuonRange_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;
class PlaneRange;

//================================================================
class HistMuCapMuonRange {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const PlaneRange& muonRange, const ClustersByPlane& globalPlaneClusters, const EventClass& evt);

  HistMuCapMuonRange()
    : doMCTruth_()
    , hMuonLastPlane_()
    , hMuonRangeGaps_()
    , hMuonMissingPlanes_()
    , hTruePProtonVsLastPlane_()
    , hTruePDeuteronVsLastPlane_()
  {}

private :
  bool doMCTruth_;

  TH1 *hMuonLastPlane_;
  TH2 *hMuonRangeGaps_;
  TH1 *hMuonMissingPlanes_;

  TH2 *hTruePProtonVsLastPlane_;
  TH2 *hTruePDeuteronVsLastPlane_;
};

#endif/*HistMuCapMuonRange_h*/
