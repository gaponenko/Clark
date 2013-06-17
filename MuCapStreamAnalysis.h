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
    ax->SetBinLabel(1+CUT_NUMAFTERTRIG, "num after trig");
    ax->SetBinLabel(1+CUT_WINTIME, "win time");
    ax->SetBinLabel(1+CUT_TGT_START, "tgt start");
    ax->SetBinLabel(1+CUT_MAX_RANGE, "MAX_RANGE");
    ax->SetBinLabel(1+CUT_RANGE_GAPS, "RANGE_GAPS");
    ax->SetBinLabel(1+CUT_OTHER_SIDE, "other side");
    ax->SetBinLabel(1+CUTS_LOOSE_PROTONS, "Loose protons");
    ax->SetBinLabel(1+CUT_MIN_RANGE, "Min range");
    ax->SetBinLabel(1+CUT_REXT, "Rext");
    ax->SetBinLabel(1+CUTS_TIGHT_PROTONS, "Tight protons");
  }

public:
  enum EventCutNumber {
    CUT_NUMAFTERTRIG,
    CUT_WINTIME,
    // UV analysis goes here
    CUT_TGT_START,
    CUT_MAX_RANGE, // in the "stream" direction
    CUT_RANGE_GAPS, // in the "stream" direction
    CUT_OTHER_SIDE,
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
    , cutMaxPlane_()
    , cutRextMax_()
    , h_cuts_r()
    , h_cuts_p()
    , hNumAfterTrigWindows_()
    , hWindowTime_()
//    , hStartPos_()
//    , hStartOffset_()
//    , hStartClusterSize_()
//    , hLastPlane_()
//    , hCCRvsPlaneDIO_()
//    , hCCRvsPlaneProtons_()
//    , hLastPlaneVsMCPstart_()
  {}

private :
  bool doMCTruth_;
  TimeWindow::StreamType cutStream_;

  double cutWinTimeMin_;
  double cutWinTimeMax_;

  int cutMaxPlane_;
  double cutRextMax_;

  TH1 *h_cuts_r;
  TH1 *h_cuts_p;

  TH1 *hNumAfterTrigWindows_;
  TH1 *hWindowTime_;

//  TH2 *hNumClusters_;
//
//  TH2 *hStartPos_;
//  TH2 *hStartOffset_; // from the muon stop position
//  TH2 *hStartClusterSize_;
//
//  TH1 *hLastPlane_;
//
//  TH2 *hCCRvsPlaneDIO_;
//  TH2 *hCCRvsPlaneProtons_;
//
//  TH2 *hLastPlaneVsMCPstart_;

  HistTDCWidth hwidthPCTightProtons_;
  HistTDCWidth hwidthDCTightProtons_;
  HistTDCWidth hwidthPCTightDIO_;
  HistTDCWidth hwidthDCTightDIO_;

  HistTDCParticleClassifier hcLooseProtons_;
  HistTDCParticleClassifier hcTightProtons_;
  HistTDCParticleClassifier hcdio_;

  MuCapUVAnalysis uvan_;
  //HistStreamAnalysis hpw_;

//  MuCapContainmentCheck rcheckProtonCandidates_;
//  HistMuCapRTruth hrtruth_;
//  HistMuCapTruth htruthLoose_;
//  HistMuCapTruth htruthPC8_;
//  HistMuCapTruth htruthMinRange_;
//  HistMuCapTruth htruthTight_;

  EventCutNumber analyze(const EventClass& evt,
                         const TimeWindowingResults& wres,
                         const ROOT::Math::XYPoint& muStopUV,
                         const std::vector<ClustersByPlane>& afterTrigGlobalClusters);
};

#endif/*MuCapStreamAnalysis_h*/
