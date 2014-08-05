// Andrei Gaponenko, 2014

#ifndef HistRangePID_h
#define HistRangePID_h

#include "WireCluster.h"

class TH1;
class TH2;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistRangePID {
public:
  void init(const std::string& hdir,
            HistogramFactory &hf,
            const ConfigFile &conf);

  void fill(const EventClass& evt,
            int itrack,
            const ClustersByPlane& protonGlobalClusters);

  HistRangePID()
    : planeRangeVsPz_(), planeRangecosVsP_()
    , trackRangeVsPz_(), trackRangecosVsP_()
  {}

private :
  TH2 *planeRangeVsPz_;
  TH2 *planeRangecosVsP_;
  TH2 *trackRangeVsPz_;
  TH2 *trackRangecosVsP_;
};

#endif/*HistRangePID_h*/
