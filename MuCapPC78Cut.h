// A cut to reduce electron contamination in the capture analysis
//
// Andrei Gaponenko, 2014

#ifndef MuCapPC78Cut_h
#define MuCapPC78Cut_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class MuCapPC78Cut {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  bool accepted(const EventClass& evt, const ClustersByPlane& globalPlaneClusters);

  bool accepted(const TDCHitWPPtr& pc7WidestHit, const TDCHitWPPtr& pc8WidestHit);

  // This one is intended for events with no hits in PC8
  bool accepted(const TDCHitWPPtr& pc7WidestHit);

private :
  double cutPC7width_;
  double cutPChalfsumWidth_;

  TH1* pc7width_before_;  // only filled for events with no pc8
  TH1* pchalfsumwidth_before_;
  TH2* pc78width_before_;

  TH1* pc7width_after_;  // only filled for events with no pc8
  TH1* pchalfsumwidth_after_;
  TH2* pc78width_after_;
};

#endif/*MuCapPC78Cut_h*/
