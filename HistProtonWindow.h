// Andrei Gaponenko, 2013

#ifndef HistProtonWindow_h
#define HistProtonWindow_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistProtonWindow {
public:
  void init(HistogramFactory &hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile &conf);

  void fill(const ClustersByPlane& globalPlaneClusters,
            const EventClass& evt);

  HistProtonWindow()
    : hClustersVsPlane_()
    , hMaxClustersVsLastPlane_()
    , hClustersVsRemainingRange_()
  {}

private :
  TH2 *hClustersVsPlane_;
  TH2 *hMaxClustersVsLastPlane_;
  TH2 *hClustersVsRemainingRange_;
};

#endif/*HistProtonWindow_h*/
