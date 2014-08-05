// Andrei Gaponenko, 2014

#ifndef HistTDCSinglePlanePID_h
#define HistTDCSinglePlanePID_h

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;


struct TDCPlanePIDResult {
  bool analyzed;
  int  nwires;
  double raw;
  double calibrated;
  TDCPlanePIDResult();
};

//================================================================
class HistTDCSinglePlanePID {
public:
  void init(const std::string& hdir,
            int globalPlaneNumber, // to use for PID
            HistogramFactory &hf,
            const ConfigFile &conf);

  TDCPlanePIDResult fill(const EventClass& evt,
                         int itrack,
                         const ClustersByPlane& protonGlobalClusters);

  HistTDCSinglePlanePID()
    : globalPlaneNumber_() , hNumClusters_(), hSingleClusterSize_()
    , hsumwcos_vs_p_1_(), hsumwcos_vs_p_2_(), hsumwcos_vs_p_3_()
  {}

private :
  int globalPlaneNumber_;
  TH1 *hNumClusters_;
  TH1 *hSingleClusterSize_;
  TH2 *hsumwcos_vs_p_1_;
  TH2 *hsumwcos_vs_p_2_;
  TH2 *hsumwcos_vs_p_3_;
};

#endif/*HistTDCSinglePlanePID_h*/
