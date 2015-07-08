// Various plots to understand the range vs p_rec discrepancy in the
// contained channel.
//
// Andrei Gaponenko, 2015

#ifndef HistTrkClass_h
#define HistTrkClass_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;
class TProfile2D;
class TEfficiency;

class HistogramFactory;
class ConfigFile;
class DetectorGeo;
class EventClass;

//================================================================
class HistTrkClass {
public:
  void init(HistogramFactory& hf,
            const std::string& hdir,
            const DetectorGeo& geom,
            const ConfigFile& conf);

  void fill(const EventClass& evt, int itrack, const ClustersByPlane& globalPlaneClusters);

private :
  TH1 *hptot_all_;

  TH1 *hptot_zcrc_;
  TH1 *hptot_zcrn_;
  TH1 *hptot_znrc_;
  TH1 *hptot_znrn_;

  TH2 *hpcos_all_;

  TH2 *hpcos_zcrc_;
  TH2 *hpcos_zcrn_;
  TH2 *hpcos_znrc_;
  TH2 *hpcos_znrn_;
};

#endif/*HistTrkClass_h*/
