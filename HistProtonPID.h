// Andrei Gaponenko, 2014

#ifndef HistProtonPID_h
#define HistProtonPID_h

#include <string>

#include "WireCluster.h"

#include "HistTDCSinglePlanePID.h"
#include "HistRangePID.h"

class TH1;
class TH2;
class TH3;
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
    , hsum78cos_vs_p_11_()
    , hsum78cos_vs_p_12_()
    , hsum78cos_vs_p_21_()
    , hsum78cos_vs_p_22_()
    , hcs78cos_vs_p_11_()
    , hcs78cos_vs_p_12_()
    , hcs78cos_vs_p_21_()
    , hcs78cos_vs_p_22_()
    , hcsPC8vsPC7vsp_()
    , hcsDC24vsDC23vsp_()
    , hcsDCAvgvsPCavgVsp_()
    , hcsDCAvgvsPCavg_p100_()
  {}

private :
  HistTDCSinglePlanePID pidPC7_;
  HistTDCSinglePlanePID pidPC8_;
  HistTDCSinglePlanePID pidDC23_;
  HistTDCSinglePlanePID pidDC24_;

  HistRangePID pidRange_;

  //----------------
  TH2 *hNumClusters78_;
  TH2 *hSingleClusterSize78_;

  TH2 *hsum78cos_vs_p_11_;
  TH2 *hsum78cos_vs_p_12_;
  TH2 *hsum78cos_vs_p_21_;
  TH2 *hsum78cos_vs_p_22_;

  // Calibrated PC width
  TH2 *hcs78cos_vs_p_11_;
  TH2 *hcs78cos_vs_p_12_;
  TH2 *hcs78cos_vs_p_21_;
  TH2 *hcs78cos_vs_p_22_;

  // DC vars
  TH2 *hNumClusters2324_;
  TH2 *hSingleClusterSize2324_;

  TH2 *hsum2324cos_vs_p_11_;
  TH2 *hsum2324cos_vs_p_12_;
  TH2 *hsum2324cos_vs_p_21_;
  TH2 *hsum2324cos_vs_p_22_;

  TH3* hcsPC8vsPC7vsp_;
  TH3* hcsDC24vsDC23vsp_;
  TH3* hcsDCAvgvsPCavgVsp_;

  TH2* hcsPC8vsPC7_p100_;
  TH2* hcsDC24vsDC23_p100_;
  TH2* hcsDCAvgvsPCavg_p100_;
};

#endif/*HistProtonPID_h*/
