// Andrei Gaponenko, 2014

#ifndef HistProtonPID_h
#define HistProtonPID_h

#include <string>

#include "WireCluster.h"

class TH1;
class TH2;
class TProfile2D;

class HistogramFactory;
class ConfigFile;
class EventClass;

//================================================================
class HistProtonPID {
public:
  void init(const std::string& hdir,
            HistogramFactory& hf,
            const ConfigFile& conf);

  void fill(const EventClass& evt, int itrack, const ClustersByPlane& protonGlobalClusters);

  HistProtonPID()
    : hNumClusters78_()
    , hSingleClusterSize78_()
    , hsum78cos_vs_p_all_()
    , hsum78cos_vs_p_11_()
    , hsum78cos_vs_p_12_()
    , hsum78cos_vs_p_21_()
    , hsum78cos_vs_p_22_()
  {}

private :
  TH2 *hNumClusters78_;
  TH2 *hSingleClusterSize78_;
  TH2 *hsum78cos_vs_p_all_;
  TH2 *hsum78cos_vs_p_11_;
  TH2 *hsum78cos_vs_p_12_;
  TH2 *hsum78cos_vs_p_21_;
  TH2 *hsum78cos_vs_p_22_;
};

#endif/*HistProtonPID_h*/
