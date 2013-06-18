// Andrei Gaponenko, 2013

#ifndef MuCapStreamAnalysis_h
#define MuCapStreamAnalysis_h

#include <string>

#include "EventClass.h"
#include "WireCluster.h"
#include "TimeWindow.h"
#include "TDCHitWP.h"

#include "HistTDCWidth.h"
//#include "HistStreamAnalysis.h"
#include "MuCapUVAnalysis.h"
#include "MuCapContainmentCheck.h"
#include "HistMuCapRTruth.h"
#include "HistMuCapTruth.h"
#include "HistTDCParticleClassifier.h"
#include "HistHitStream.h"

#include "Math/Point2D.h"
#include "TAxis.h"
class TH1;
class TH2;

class HistogramFactory;
class DetectorGeo;
class ConfigFile;

//================================================================
class MuCapStreamAnalysis {
  void set_cut_bin_labels(TAxis* ax) {
    ax->SetBinLabel(1+CUT_NOHITS, "No hits");
    ax->SetBinLabel(1+CUT_MULTIWIN, "Multiple windows");
    ax->SetBinLabel(1+CUT_WINTIME, "win time");
    ax->SetBinLabel(1+CUT_Z_CONTAINED, "Z contained");
    //ax->SetBinLabel(1+CUT_TGT_START, "tgt start");
    //ax->SetBinLabel(1+CUT_MAX_RANGE, "MAX_RANGE");
    //ax->SetBinLabel(1+CUT_RANGE_GAPS, "RANGE_GAPS");
    //ax->SetBinLabel(1+CUT_OTHER_SIDE, "other side");
    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");
    ax->SetBinLabel(1+CUT_MIN_RANGE, "Min range");
    ax->SetBinLabel(1+CUT_REXT, "Rext");
    ax->SetBinLabel(1+CUTS_TIGHT_PROTONS, "Tight protons");
  }

public:
  enum EventCutNumber {
    CUT_NOHITS,
    CUT_MULTIWIN,
    CUT_WINTIME,
    CUT_Z_CONTAINED,
    // UV analysis goes here
    CUT_TGT_START,
    CUT_RANGE_GAPS, // in the "stream" direction
    CUT_OTHER_SIDE,

    CUT_MAX_RANGE, // in the "stream" direction

    CUTS_LOOSE_PROTONS,
    CUT_MIN_RANGE,
    CUT_REXT,
    CUTS_TIGHT_PROTONS,
    CUTS_END
  };

  void init(HistogramFactory &hf, const std::string& hdir,
            const DetectorGeo& geom, const ConfigFile &conf,
            TimeWindow::StreamType cutWinStream, double cutWinTimeMin);

  void process(const EventClass& evt,
               const TimeWindowingResults& wres,
               const ROOT::Math::XYPoint& muStopUV,
               const std::vector<ClustersByPlane>& afterTrigGlobalClusters);

  MuCapStreamAnalysis()
    : doMCTruth_(false)
    , cutStream_(TimeWindow::DOWNSTREAM)
    , cutWinTimeMin_()
    , cutWinTimeMax_()
    , cutZContainedNumPlanes_()
    , cutRextMax_()
    , h_cuts_r()
    , h_cuts_p()
    , hMultiWindowHits_()
    , hMultiWindowRanges_()
    , hNumAfterTrigWindows_()
    , hWindowTimeBefore_()
    , hWindowTimeAfter_()
    , hNumVetoHits_()
  {}

private :
  bool doMCTruth_;
  TimeWindow::StreamType cutStream_;

  double cutWinTimeMin_;
  double cutWinTimeMax_;

  // The minimal required number of hit-free planes
  // on both ends to treat an event as "z contained"
  int cutZContainedNumPlanes_;

  double cutRextMax_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH2 *hMultiWindowHits_;
  TH2 *hMultiWindowRanges_;
  TH1 *hNumAfterTrigWindows_;

  TH1 *hWindowTimeBefore_;
  TH1 *hWindowTimeAfter_;
  HistHitStream hhsAfterTimeCuts_;

  TH1 *hNumVetoHits_;
  HistHitStream hhsZContained_;

  HistTDCWidth hwidthPCTightProtons_;
  HistTDCWidth hwidthDCTightProtons_;
  HistTDCWidth hwidthPCTightDIO_;
  HistTDCWidth hwidthDCTightDIO_;

  HistTDCParticleClassifier hcLooseProtons_;
  HistTDCParticleClassifier hcTightProtons_;
  HistTDCParticleClassifier hcdio_;

  MuCapUVAnalysis uvan_;

  EventCutNumber analyze(const EventClass& evt,
                         const TimeWindowingResults& wres,
                         const ROOT::Math::XYPoint& muStopUV,
                         const std::vector<ClustersByPlane>& afterTrigGlobalClusters);
};

#endif/*MuCapStreamAnalysis_h*/
