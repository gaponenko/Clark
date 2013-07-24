// Andrei Gaponenko, 2013

#ifndef HistMuCapFinal_h
#define HistMuCapFinal_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistMuCapFinal {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const ClustersByPlane& globalPlaneClusters,
            const EventClass& evt);

  HistMuCapFinal()
    : hLastPlane_()
    , hNumPlanesVsWires_()
  {}

private :
  TH1 *hLastPlane_;
  TH2 *hNumPlanesVsWires_;
};

#endif/*HistMuCapFinal_h*/
